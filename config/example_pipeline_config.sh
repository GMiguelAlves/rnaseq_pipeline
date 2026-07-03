#!/usr/bin/env bash

# Backward-compatible example file.
# Prefer copying config/user_settings_template.sh to config/user_settings.sh.

# shellcheck source=user_settings_template.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/user_settings_template.sh"
