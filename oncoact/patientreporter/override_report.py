import argparse
import json
from google.cloud import pubsub
from cli_util import perform_prod_test

TOPIC = 'foo'


def main():
    """
    This script is used to easily emit an override report command/event to the reporting pipeline.
    """
    parser = argparse.ArgumentParser(description='Emit a report override event to the reporting pipeline')
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    perform_prod_test(args.profile)

    project = 'hmf-pipeline-development'  # TODO set dynamically
    assemble_and_emit_event(args.tumor_sample_barcode, project)


def assemble_and_emit_event(tumor_sample_barcode, project):
    """
    Gather all the input information needed to emit the report override event.

    :param tumor_sample_barcode: the tumor sample barcode of the report to emit the event for.
    :param project: the gcp project to emit this event to.
    """
    special_remark = input('Enter a special remark...\n')
    remark_is_external = input('Is remark external? y/n\n').lower() == 'y'
    rose_override = input('Enter rose override...\n')
    comments = input('Enter comments...\n')

    client = pubsub.PublisherClient()
    topic_path = client.topic_path(project=project, topic='report.override')

    res = emit_correction_event(client,
                                topic_path,
                                tumor_sample_barcode,
                                special_remark,
                                remark_is_external,
                                rose_override,
                                comments)

    print(f'Event emitted! (message id: {res})')
    exit(0)


def emit_correction_event(client: pubsub.PublisherClient,
                          topic_path: str,
                          tumor_sample_barcode: str,
                          special_remark: str,
                          remark_is_external: bool,
                          rose_override: str,
                          comments: str) -> str:
    """
    Creates and emits the report override event.

    :param client: the pubsub client.
    :param topic_path: the topic path (which includes both project and topic).
    :param tumor_sample_barcode: the tumor sample barcode to emit the event for.
    :param special_remark: the special remark section.
    :param remark_is_external: boolean value that indicates if this remark is external.
    :param rose_override: the rose override section.
    :param comments: the comments section.
    :return: the id of the emitted event.
    """
    data = {
        'tumorSampleBarcode': tumor_sample_barcode,
        'specialRemark': special_remark,
        'remarkIsExternal': remark_is_external,
        'roseOverride': rose_override,
        'comments': comments
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode(encoding='utf-8')

    return client.publish(topic=topic_path, data=encoded_message).result()


if __name__ == '__main__':
    main()
