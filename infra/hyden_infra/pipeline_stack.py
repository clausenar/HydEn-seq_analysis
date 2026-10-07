"""HydEn-seq pipeline infrastructure as code.

Defines the same event-driven pipeline that was first provisioned by hand
(see the "AWS deployment" section of the top-level README):

    S3 landing/{sample}/*.fastq  --(ObjectCreated)-->  Lambda landing watcher
        --(DynamoDB conditional claim)-->  batch:SubmitJob
        -->  AWS Batch on Fargate running the pipeline image from ECR
        -->  output/{sample}/ in S3, raw input archived to processed_raw/{sample}/

Each stage ("dev", "prod", ...) gets its own copy of every resource, so a
change can be deployed and tested next to the live pipeline before it
replaces it. The container image is not built here: GitHub Actions builds it
and pushes it to the existing ECR repository, which this stack only imports.
"""

from pathlib import Path

from aws_cdk import (
    Duration,
    RemovalPolicy,
    Size,
    Stack,
    aws_batch as batch,
    aws_cloudwatch as cloudwatch,
    aws_cloudwatch_actions as cw_actions,
    aws_dynamodb as dynamodb,
    aws_ec2 as ec2,
    aws_ecr as ecr,
    aws_ecs as ecs,
    aws_events as events,
    aws_events_targets as targets,
    aws_iam as iam,
    aws_lambda as lambda_,
    aws_logs as logs,
    aws_s3 as s3,
    aws_s3_notifications as s3n,
    aws_sns as sns,
    aws_sns_subscriptions as subs,
)
from constructs import Construct

LAMBDA_SOURCE = Path(__file__).resolve().parents[2] / "lambda"

LANDING_PREFIX = "landing/"
PROCESSED_RAW_PREFIX = "processed_raw/"


class HydenSeqPipelineStack(Stack):
    def __init__(
        self,
        scope: Construct,
        construct_id: str,
        *,
        stage: str,
        ecr_repository_name: str = "hyden-seq-pipeline",
        image_tag: str = "latest",
        alert_email: str | None = None,
        **kwargs,
    ) -> None:
        super().__init__(scope, construct_id, **kwargs)

        # Production keeps its data if the stack is deleted; any other stage
        # is disposable, so `cdk destroy` leaves nothing behind (and nothing billed).
        is_prod = stage == "prod"
        removal = RemovalPolicy.RETAIN if is_prod else RemovalPolicy.DESTROY

        # ---------------------------------------------------------------- storage
        self.bucket = s3.Bucket(
            self,
            "DataBucket",
            block_public_access=s3.BlockPublicAccess.BLOCK_ALL,
            encryption=s3.BucketEncryption.S3_MANAGED,
            enforce_ssl=True,
            removal_policy=removal,
            auto_delete_objects=not is_prod,
            lifecycle_rules=[
                # Raw fastqs are only re-read for reprocessing, so move them to
                # cheaper storage once they have been processed for a month.
                s3.LifecycleRule(
                    id="ArchiveProcessedRaw",
                    prefix=PROCESSED_RAW_PREFIX,
                    transitions=[
                        s3.Transition(
                            storage_class=s3.StorageClass.INFREQUENT_ACCESS,
                            transition_after=Duration.days(30),
                        )
                    ],
                ),
                s3.LifecycleRule(
                    id="AbortIncompleteUploads",
                    abort_incomplete_multipart_upload_after=Duration.days(7),
                ),
            ],
        )

        # Exactly-once job submission + job status per sample (see the Lambda docstring).
        self.status_table = dynamodb.Table(
            self,
            "SampleStatusTable",
            partition_key=dynamodb.Attribute(name="sample_id", type=dynamodb.AttributeType.STRING),
            billing_mode=dynamodb.BillingMode.PAY_PER_REQUEST,
            point_in_time_recovery_specification=dynamodb.PointInTimeRecoverySpecification(
                point_in_time_recovery_enabled=is_prod
            ),
            removal_policy=removal,
        )

        # ---------------------------------------------------------------- network
        # Public subnets only and no NAT gateway: Fargate tasks get a public IP
        # to reach ECR, which avoids the ~$30/month a NAT gateway would cost.
        # S3 traffic goes through a free gateway endpoint instead of the internet.
        vpc = ec2.Vpc(
            self,
            "Vpc",
            # The account's AZ list is pinned in cdk.json context, so `cdk synth`
            # needs no AWS credentials (no lookup).
            max_azs=2,
            nat_gateways=0,
            subnet_configuration=[
                ec2.SubnetConfiguration(name="public", subnet_type=ec2.SubnetType.PUBLIC)
            ],
            gateway_endpoints={"S3": ec2.GatewayVpcEndpointOptions(service=ec2.GatewayVpcEndpointAwsService.S3)},
        )

        # ---------------------------------------------------------------- compute
        compute_env = batch.FargateComputeEnvironment(
            self,
            "FargateComputeEnv",
            vpc=vpc,
            vpc_subnets=ec2.SubnetSelection(subnet_type=ec2.SubnetType.PUBLIC),
            maxv_cpus=16,
        )

        self.job_queue = batch.JobQueue(
            self,
            "JobQueue",
            compute_environments=[
                batch.OrderedComputeEnvironment(compute_environment=compute_env, order=1)
            ],
        )

        # The container's own role: read and write the data bucket only
        # (write + delete are needed for `aws s3 mv` to processed_raw/).
        job_role = iam.Role(
            self,
            "BatchJobRole",
            assumed_by=iam.ServicePrincipal("ecs-tasks.amazonaws.com"),
            description="HydEn-seq pipeline container: read/write its own data bucket only",
        )
        self.bucket.grant_read_write(job_role)

        repository = ecr.Repository.from_repository_name(self, "PipelineRepo", ecr_repository_name)

        job_logs = logs.LogGroup(
            self,
            "BatchJobLogs",
            retention=logs.RetentionDays.ONE_MONTH,
            removal_policy=removal,
        )

        self.job_definition = batch.EcsJobDefinition(
            self,
            "PipelineJobDefinition",
            container=batch.EcsFargateContainerDefinition(
                self,
                "PipelineContainer",
                image=ecs.ContainerImage.from_ecr_repository(repository, image_tag),
                cpu=8,
                memory=Size.gibibytes(16),
                assign_public_ip=True,
                job_role=job_role,
                logging=ecs.LogDriver.aws_logs(stream_prefix="hyden-seq", log_group=job_logs),
            ),
            retry_attempts=1,
            timeout=Duration.hours(12),
        )

        # ------------------------------------------------------- landing watcher
        watcher_logs = logs.LogGroup(
            self,
            "LandingWatcherLogs",
            retention=logs.RetentionDays.ONE_MONTH,
            removal_policy=removal,
        )

        self.landing_watcher = lambda_.Function(
            self,
            "LandingWatcher",
            runtime=lambda_.Runtime.PYTHON_3_12,
            handler="landing_watcher_lambda.handler",
            code=lambda_.Code.from_asset(
                str(LAMBDA_SOURCE),
                exclude=["test_*", "requirements-dev.txt", "__pycache__", "*.zip"],
            ),
            timeout=Duration.seconds(30),
            memory_size=256,
            log_group=watcher_logs,
            environment={
                "STATUS_TABLE_NAME": self.status_table.table_name,
                "BATCH_JOB_QUEUE": self.job_queue.job_queue_arn,
                "BATCH_JOB_DEFINITION": self.job_definition.job_definition_arn,
            },
        )

        # Least privilege, matching what the Lambda actually calls:
        # s3:ListBucket/GetObject, three DynamoDB item actions, batch:SubmitJob.
        self.bucket.grant_read(self.landing_watcher)
        self.status_table.grant(
            self.landing_watcher,
            "dynamodb:GetItem",
            "dynamodb:PutItem",
            "dynamodb:UpdateItem",
        )
        self.landing_watcher.add_to_role_policy(
            iam.PolicyStatement(
                actions=["batch:SubmitJob"],
                resources=[self.job_queue.job_queue_arn, self.job_definition.job_definition_arn],
            )
        )

        self.bucket.add_event_notification(
            s3.EventType.OBJECT_CREATED,
            s3n.LambdaDestination(self.landing_watcher),
            s3.NotificationKeyFilter(prefix=LANDING_PREFIX),
        )

        # ---------------------------------------------------------------- alerting
        self.alerts = sns.Topic(self, "PipelineAlerts", display_name=f"HydEn-seq {stage} alerts")
        if alert_email:
            self.alerts.add_subscription(subs.EmailSubscription(alert_email))

        # Batch has no built-in "failed jobs" metric, so failed jobs on this
        # queue are routed to the alert topic through EventBridge.
        events.Rule(
            self,
            "FailedJobRule",
            description="Alert when a HydEn-seq Batch job fails",
            event_pattern=events.EventPattern(
                source=["aws.batch"],
                detail_type=["Batch Job State Change"],
                detail={"status": ["FAILED"], "jobQueue": [self.job_queue.job_queue_arn]},
            ),
            targets=[targets.SnsTopic(self.alerts)],
        )

        watcher_errors = cloudwatch.Alarm(
            self,
            "LandingWatcherErrors",
            alarm_description="The landing watcher Lambda raised an error",
            metric=self.landing_watcher.metric_errors(period=Duration.minutes(5)),
            threshold=1,
            evaluation_periods=1,
            treat_missing_data=cloudwatch.TreatMissingData.NOT_BREACHING,
        )
        watcher_errors.add_alarm_action(cw_actions.SnsAction(self.alerts))
