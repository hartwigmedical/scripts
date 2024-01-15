import json
import argparse
from google.cloud import pubsub

from input_util import multi_line_input


def main():
    """
    This script is used to easily emit a fail report command/event to the reporting pipeline.
    """
    parser = argparse.ArgumentParser(description="Emit a QC fail event to the reporting pipeline.")
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'preview', 'prod'], default='pilot')
    args = parser.parse_args()

    tumor_sample_barcode = args.tumor_sample_barcode

    if args.profile == 'pilot':
        project = 'hmf-pipeline-development'
    elif args.profile == 'prod' or args.profile == 'preview':
        project = 'hmf-pipeline-prod-e45b00f2'
    else:
        raise ValueError(f"Profile '{args.profile}' not recognized.")
    assemble_and_emit_qc_fail_event(tumor_sample_barcode, project)


def assemble_and_emit_qc_fail_event(tumor_sample_barcode: str,
                                    project: str):
    reason = _prompt_user_for_reason()
    fail_reason_comment = multi_line_input("Please add an optional fail reason comment\n")
    fail_reason_comment = fail_reason_comment if fail_reason_comment else None

    add_correction = input('Add correction? y/n\n').lower() == 'y'
    correction = None
    if add_correction:
        remark_is_external = input('Is remark external? y/n\n').lower() == 'y'
        comments = multi_line_input('Enter comments...\n')
        correction = {
            'comments': comments,
            'remarkIsExternal': remark_is_external
        }

    data = {
        'tumorSampleBarcode': tumor_sample_barcode,
        'reason': reason,
        'failReasonComment': fail_reason_comment,
        'correction': correction
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode('utf-8')

    pubsub_client = pubsub.PublisherClient()
    topic_path = pubsub_client.topic_path(project=project, topic='report.fail')

    message_id = pubsub_client.publish(topic=topic_path, data=encoded_message).result()
    print(f"event emitted! (message id: '{message_id}')")


def _prompt_user_for_reason():
    choices = ['wgs_processing_issue', 'wgs_isolation_fail', 'wgs_tcp_shallow_fail', 'wgs_preparation_fail',
               'wgs_tumor_processing_issue', 'wgs_pipeline_fail', 'wgs_tcp_fail', 'panel_result_report_fail']
    print('Please provide a reason for the failure, the possible reasons are:')
    print('\n'.join(choices))
    reason = input()
    if reason not in choices:
        print("not a valid option\n")
        return _prompt_user_for_reason()
    return reason.upper()


if __name__ == '__main__':
    main()
