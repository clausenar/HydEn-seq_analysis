"""Tests for the synthesised CloudFormation template.

These run `cdk synth` in-process and check the template, so they need no AWS
credentials and touch no real account. Run from infra/:

    .venv/bin/pytest -v
"""

import json
from pathlib import Path

import aws_cdk as cdk
import pytest
from aws_cdk.assertions import Match, Template

from hyden_infra.pipeline_stack import HydenSeqPipelineStack

CDK_CONTEXT = json.loads((Path(__file__).resolve().parents[1] / "cdk.json").read_text())["context"]
ENV = cdk.Environment(account="069509443906", region="us-east-1")


def synth(stage="dev", **kwargs):
    app = cdk.App(context=CDK_CONTEXT)
    stack = HydenSeqPipelineStack(app, f"Test-{stage}", stage=stage, env=ENV, **kwargs)
    return Template.from_stack(stack)


@pytest.fixture(scope="module")
def dev():
    return synth("dev")


@pytest.fixture(scope="module")
def prod():
    return synth("prod")


def policy_actions(template, role_logical_prefix):
    """All IAM actions granted to the role whose logical ID starts with the prefix."""
    actions = set()
    for logical_id, policy in template.find_resources("AWS::IAM::Policy").items():
        if not logical_id.startswith(role_logical_prefix):
            continue
        for statement in policy["Properties"]["PolicyDocument"]["Statement"]:
            action = statement["Action"]
            actions.update(action if isinstance(action, list) else [action])
    return actions


# ------------------------------------------------------------------ storage
def test_bucket_is_private_encrypted_and_https_only(dev):
    dev.has_resource_properties(
        "AWS::S3::Bucket",
        {
            "PublicAccessBlockConfiguration": {
                "BlockPublicAcls": True,
                "BlockPublicPolicy": True,
                "IgnorePublicAcls": True,
                "RestrictPublicBuckets": True,
            },
            "BucketEncryption": Match.object_like({}),
        },
    )
    dev.has_resource_properties(
        "AWS::S3::BucketPolicy",
        {
            "PolicyDocument": {
                "Statement": Match.array_with(
                    [
                        Match.object_like(
                            {"Effect": "Deny", "Condition": {"Bool": {"aws:SecureTransport": "false"}}}
                        )
                    ]
                )
            }
        },
    )


def test_processed_raw_moves_to_cheaper_storage(dev):
    dev.has_resource_properties(
        "AWS::S3::Bucket",
        {
            "LifecycleConfiguration": {
                "Rules": Match.array_with(
                    [
                        Match.object_like(
                            {
                                "Prefix": "processed_raw/",
                                "Transitions": [{"StorageClass": "STANDARD_IA", "TransitionInDays": 30}],
                            }
                        )
                    ]
                )
            }
        },
    )


def test_status_table_matches_lambda_contract(dev):
    dev.has_resource_properties(
        "AWS::DynamoDB::Table",
        {
            "KeySchema": [{"AttributeName": "sample_id", "KeyType": "HASH"}],
            "BillingMode": "PAY_PER_REQUEST",
        },
    )


def test_prod_keeps_data_and_dev_is_disposable(dev, prod):
    prod.has_resource("AWS::S3::Bucket", {"DeletionPolicy": "Retain"})
    prod.has_resource("AWS::DynamoDB::Table", {"DeletionPolicy": "Retain"})
    prod.resource_count_is("Custom::S3AutoDeleteObjects", 0)

    dev.has_resource("AWS::S3::Bucket", {"DeletionPolicy": "Delete"})
    dev.has_resource("AWS::DynamoDB::Table", {"DeletionPolicy": "Delete"})


# ------------------------------------------------------------------ compute
def test_job_definition_runs_on_fargate_with_8_vcpu_16_gb(dev):
    dev.has_resource_properties(
        "AWS::Batch::JobDefinition",
        {
            "PlatformCapabilities": ["FARGATE"],
            "ContainerProperties": Match.object_like(
                {
                    "ResourceRequirements": Match.array_with(
                        [{"Type": "MEMORY", "Value": "16384"}, {"Type": "VCPU", "Value": "8"}]
                    )
                }
            ),
        },
    )


def test_no_nat_gateway(dev):
    # A NAT gateway would cost ~$30/month even while the pipeline is idle.
    dev.resource_count_is("AWS::EC2::NatGateway", 0)


def test_compute_scales_to_zero_with_a_cap(dev):
    dev.has_resource_properties(
        "AWS::Batch::ComputeEnvironment",
        {"ComputeResources": Match.object_like({"Type": "FARGATE", "MaxvCpus": 16})},
    )


# ---------------------------------------------------------- landing watcher
def test_landing_watcher_gets_its_configuration(dev):
    dev.has_resource_properties(
        "AWS::Lambda::Function",
        {
            "Handler": "landing_watcher_lambda.handler",
            "Environment": {
                "Variables": {
                    "STATUS_TABLE_NAME": Match.any_value(),
                    "BATCH_JOB_QUEUE": Match.any_value(),
                    "BATCH_JOB_DEFINITION": Match.any_value(),
                }
            },
        },
    )


def test_only_landing_uploads_trigger_the_watcher(dev):
    dev.has_resource_properties(
        "Custom::S3BucketNotifications",
        {
            "NotificationConfiguration": {
                "LambdaFunctionConfigurations": [
                    Match.object_like(
                        {
                            "Events": ["s3:ObjectCreated:*"],
                            "Filter": {"Key": {"FilterRules": [{"Name": "prefix", "Value": "landing/"}]}},
                        }
                    )
                ]
            }
        },
    )


def test_landing_watcher_is_least_privilege(dev):
    actions = policy_actions(dev, "LandingWatcherServiceRole")
    assert actions == {
        "batch:SubmitJob",
        "dynamodb:GetItem",
        "dynamodb:PutItem",
        "dynamodb:UpdateItem",
        "s3:GetBucket*",
        "s3:GetObject*",
        "s3:List*",
    }
    # It must never write or delete data, or touch other Batch jobs.
    assert not any(a.startswith(("s3:Put", "s3:Delete")) for a in actions)
    assert "batch:*" not in actions and "dynamodb:*" not in actions


def test_container_can_read_write_its_bucket_only(dev):
    actions = policy_actions(dev, "BatchJobRole")
    assert {"s3:GetObject*", "s3:PutObject", "s3:DeleteObject*"} <= actions
    assert not any(a.split(":")[0] not in {"s3"} for a in actions)


# ----------------------------------------------------------------- alerting
def test_failed_jobs_are_alerted(dev):
    dev.has_resource_properties(
        "AWS::Events::Rule",
        {
            "EventPattern": Match.object_like(
                {
                    "source": ["aws.batch"],
                    "detail-type": ["Batch Job State Change"],
                    "detail": Match.object_like({"status": ["FAILED"]}),
                }
            ),
            "Targets": [Match.object_like({"Arn": Match.any_value()})],
        },
    )
    dev.resource_count_is("AWS::CloudWatch::Alarm", 1)


def test_alert_email_is_optional():
    synth("dev").resource_count_is("AWS::SNS::Subscription", 0)
    synth("dev", alert_email="someone@example.com").resource_count_is("AWS::SNS::Subscription", 1)
