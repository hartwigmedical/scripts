import argparse
import json
from google.cloud import pubsub

TOPIC = ''


def main(tumor_isolation_barcode):
    special_remark = input('Enter a special remark...')
    remark_is_external = input('Is remark external? y/n') == 'y' or 'Y'
    rose_override = input('Enter rose override...')
    comments = input('Enter comments...')

    emit_correction_event(pubsub.PublisherClient(),
                          tumor_isolation_barcode,
                          special_remark,
                          remark_is_external,
                          rose_override,
                          comments)

    print('Event emitted!')
    exit(0)


def emit_correction_event(client: pubsub.PublisherClient,
                          tumor_isolation_barcode: str,
                          special_remark: str,
                          remark_is_external: bool,
                          rose_override: str,
                          comments: str) -> None:
    data = {
        'ReportOverride': {
            'tumorIsolationBarcode': tumor_isolation_barcode,
            'specialRemark': special_remark,
            'remarkIsExternal': remark_is_external,
            'roseOverride': rose_override,
            'comments': comments
        }
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode(encoding='utf-8')

    client.publish(topic=TOPIC, data=encoded_message, subject='report', event='override')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Emit a override event for a given report')
    parser.add_argument('tumor_isolation_barcode')
    args = parser.parse_args()
    main(args.tumor_isolation_barcode)
