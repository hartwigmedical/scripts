import json
import argparse
from google.cloud import pubsub


def main():
    parser = argparse.ArgumentParser(description="Emit a QC fail event for a given report.")
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()
    if args.profile == 'prod':
        prod_warn = input("Warning: you are running in prod. Type 'y' to continue.")
        if prod_warn.lower() != 'y':
            print('Program aborted')
            exit(1)
    tumor_sample_barcode = args.tumor_sample_barcode
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
    client = pubsub.PublisherClient()
    project = 'hmf-pipeline-development'  # TODO make dynamic
    topic_path = client.topic_path(project, 'report.fail')
    res = emit_fail_report(client, topic_path, reason, tumor_sample_barcode)
    print(f"event emitted! (message id: '{res}')")
    exit(0)


def emit_fail_report(client: pubsub.PublisherClient,
                     topic_path: str,
                     reason: str,
                     tumor_sample_barcode: str):
    data = {
        'reason': reason,
        'tumorSampleBarcode': tumor_sample_barcode
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode('utf-8')
    return client.publish(topic=topic_path, data=encoded_message).result()


if __name__ == '__main__':
    main()
