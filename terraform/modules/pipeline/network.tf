# Public subnets only and no NAT gateway: Fargate tasks get a public IP to
# reach ECR, which avoids the ~$30/month a NAT gateway would cost. S3 traffic
# goes through a free gateway endpoint instead of the internet.

resource "aws_vpc" "this" {
  cidr_block           = var.vpc_cidr
  enable_dns_support   = true
  enable_dns_hostnames = true

  tags = { Name = local.name }
}

resource "aws_internet_gateway" "this" {
  vpc_id = aws_vpc.this.id

  tags = { Name = local.name }
}

resource "aws_subnet" "public" {
  count = length(var.availability_zones)

  vpc_id            = aws_vpc.this.id
  availability_zone = var.availability_zones[count.index]
  cidr_block        = cidrsubnet(var.vpc_cidr, 8, count.index)

  tags = { Name = "${local.name}-public-${var.availability_zones[count.index]}" }
}

resource "aws_route_table" "public" {
  vpc_id = aws_vpc.this.id

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.this.id
  }

  tags = { Name = "${local.name}-public" }
}

resource "aws_route_table_association" "public" {
  count = length(aws_subnet.public)

  subnet_id      = aws_subnet.public[count.index].id
  route_table_id = aws_route_table.public.id
}

resource "aws_vpc_endpoint" "s3" {
  vpc_id            = aws_vpc.this.id
  service_name      = "com.amazonaws.${data.aws_region.current.region}.s3"
  vpc_endpoint_type = "Gateway"
  route_table_ids   = [aws_route_table.public.id]

  tags = { Name = "${local.name}-s3" }
}

# Batch tasks only make outbound calls (ECR, S3, CloudWatch Logs); nothing
# may connect in.
resource "aws_security_group" "batch" {
  name        = "${local.name}-batch"
  description = "HydEn-seq Batch tasks: outbound only"
  vpc_id      = aws_vpc.this.id

  egress {
    description = "HTTPS to ECR, S3 and CloudWatch Logs"
    from_port   = 443
    to_port     = 443
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

# Takes over the VPC's default security group and leaves it with no rules,
# so anything accidentally launched without a group gets no network access.
resource "aws_default_security_group" "this" {
  vpc_id = aws_vpc.this.id
}
