#!/bin/bash

prefix=`dirname $(readlink $0 || echo $0)`

${prefix}/pilot_do_run_base_patient_reporter \
    -not_analysable \
    -not_analysable_reason shallow_seq \
    -not_analysable_sample "$@"