import argparse
import json
from google.cloud import pubsub

TOPIC = 'foo'


def main():
    parser = argparse.ArgumentParser(description='Emit a override event for a given report')
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    if args.profile == 'prod':
        prod_warn = input("Warning: you are running in prod. Type 'y' to continue.")
        if prod_warn.lower() != 'y':
            print('Program aborted')
            exit(1)

    project = 'hmf-pipeline-development'  # TODO set dynamically
    assemble_event(args.tumor_sample_barcode, project)


def assemble_event(tumor_sample_barcode, project):
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
