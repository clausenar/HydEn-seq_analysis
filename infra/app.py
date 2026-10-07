#!/usr/bin/env python3
"""CDK entry point.

    cdk synth                       # dev stage (default)
    cdk synth -c stage=prod
    cdk deploy -c stage=dev -c alert_email=you@example.com
"""

import aws_cdk as cdk

from hyden_infra.pipeline_stack import HydenSeqPipelineStack

app = cdk.App()

stage = app.node.try_get_context("stage") or "dev"

HydenSeqPipelineStack(
    app,
    f"HydenSeqPipeline-{stage}",
    stage=stage,
    alert_email=app.node.try_get_context("alert_email"),
    env=cdk.Environment(account="069509443906", region="us-east-1"),
)

cdk.Tags.of(app).add("project", "hyden-seq")
cdk.Tags.of(app).add("stage", stage)

app.synth()
