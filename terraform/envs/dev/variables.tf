variable "alert_email" {
  description = "Optional email for pipeline alerts. Pass with -var or TF_VAR_alert_email; not committed."
  type        = string
  default     = null
}
