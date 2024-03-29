import argparse
import json
from google.cloud import pubsub

from input_util import multi_line_input


def main():
    """
    This script is used to easily emit an override report command/event to the reporting pipeline.
    """
    parser = argparse.ArgumentParser(description='Emit a report override event to the reporting pipeline')
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'preview', 'prod'], default='pilot')
    args = parser.parse_args()

    if args.profile == 'pilot':
        project = 'hmf-pipeline-development'
    elif args.profile == 'prod' or args.profile == 'preview':
        project = 'hmf-pipeline-prod-e45b00f2'
    else:
        raise ValueError(f"Profile '{args.profile}' not recognized.")
    assemble_and_emit_correction_event(args.tumor_sample_barcode, project)


def assemble_and_emit_correction_event(tumor_sample_barcode: str, project: str):
    special_remark = multi_line_input('Enter an optional special remark...\n')
    special_remark = special_remark if special_remark else None
    rose_override = multi_line_input('Enter an optional rose override...\n')
    rose_override = rose_override if rose_override else None

    correction = None
    add_correction = input('Add correction? y/n\n').lower() == 'y'
    if add_correction:
        remark_is_external = input('Is remark external? y/n\n').lower() == 'y'
        comments = multi_line_input('Enter comments...\n')
        correction = {
            'comments': comments,
            'remarkIsExternal': remark_is_external
        }

    data = {
        'tumorSampleBarcode': tumor_sample_barcode,
        'specialRemark': special_remark,
        'correction': correction,
        'roseOverride': rose_override,
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode(encoding='utf-8')
    pubsub_client = pubsub.PublisherClient()
    topic_path = pubsub_client.topic_path(project, 'report.override')

    message_id = pubsub_client.publish(topic=topic_path, data=encoded_message).result()

    print(f'Event emitted! (message id: {message_id})')


if __name__ == '__main__':
    main()
