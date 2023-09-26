import json
import argparse
from google.cloud import pubsub

TOPIC = ''


def main():
    parser = argparse.ArgumentParser(description="Emit a QC fail event for a given report.")
    parser.add_argument('run_id')
    parser.add_argument('--topic', default=TOPIC)
    args = parser.parse_args()
    run_id = int(args.run_id)
    topic = args.topic

    choices = ['wgs_processing_issue', 'wgs_isolation_fail', 'wgs_tcp_shallow_fail', 'wgs_preparation_fail',
               'wgs_tumor_processing_issue', 'wgs_pipeline_fail', 'wgs_tcp_fail']

    print('Please provide a reason for the failure, the possible reasons are:')
    print('\n'.join(choices))
    reason = input()
    if reason not in choices:
        cont = input(
            f"'{reason}' is not recognized as a valid reason, are you sure you want to proceed?"
            f" This might cause an exception (y/n)")
        if cont.lower() != 'y':
            exit(1)

    emit_fail_report(pubsub.PublisherClient(), topic, reason, run_id)
    print('event emitted!')
    exit(0)


async def emit_fail_report(client: pubsub.PublisherClient,
                           topic: str,
                           reason: str,
                           run_id: int):
    data = {
        'reason': reason,
        'runId': run_id
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode('utf-8')
    await client.publish(topic=topic, data=encoded_message, subject='report', event='fail')


if __name__ == '__main__':
    main()
