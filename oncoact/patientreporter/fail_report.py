import json
import argparse
from google.cloud import pubsub

TOPIC = ''


def main(run_id):
    reason = input('Please provide a reason for the failure')
    emit_fail_report(pubsub.PublisherClient(), reason, run_id)
    print('event emitted!')
    exit(0)


def emit_fail_report(client: pubsub.PublisherClient,
                     reason: str,
                     run_id: int):
    data = {
        'reason': reason,
        'runId': run_id
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode('utf-8')

    client.publish(topic=TOPIC, data=encoded_message, subject='report', event='fail')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Emit a QC fail event for a given report.")
    parser.add_argument('run_id')
    args = parser.parse_args()
    main(int(args.run_id))
