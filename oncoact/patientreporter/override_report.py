import argparse
import json
from google.cloud import pubsub

TOPIC = 'foo'


def main():
    """
    This script is used to easily emit an override report command/event to the reporting pipeline.
    """
    parser = argparse.ArgumentParser(description='Emit a report override event to the reporting pipeline')
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    project = 'hmf-pipeline-development'  # TODO set dynamically
    assemble_and_emit_correction_event(args.tumor_sample_barcode, project)


def assemble_and_emit_correction_event(tumor_sample_barcode: str, project: str):
    special_remark = input('Enter a special remark...\n')
    remark_is_external = input('Is remark external? y/n\n').lower() == 'y'
    rose_override = input('Enter rose override...\n')
    comments = input('Enter comments...\n')

    data = {
        'tumorSampleBarcode': tumor_sample_barcode,
        'specialRemark': special_remark,
        'remarkIsExternal': remark_is_external,
        'roseOverride': rose_override,
        'comments': comments
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode(encoding='utf-8')
    pubsub_client = pubsub.PublisherClient()
    topic_path = pubsub_client.topic_path(project, 'report.override')

    message_id = pubsub_client.publish(topic=topic_path, data=encoded_message).result()

    print(f'Event emitted! (message id: {message_id})')


if __name__ == '__main__':
    main()
