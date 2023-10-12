import json
import argparse
from google.cloud import pubsub


def main():
    """
    This script is used to easily emit a fail report command/event to the reporting pipeline.
    """
    parser = argparse.ArgumentParser(description="Emit a QC fail event to the reporting pipeline.")
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    tumor_sample_barcode = args.tumor_sample_barcode
    project = 'hmf-pipeline-development'
    assemble_and_emit_fail_report(tumor_sample_barcode, project)


def assemble_and_emit_fail_report(tumor_sample_barcode: str,
                                  project: str):
    reason = _prompt_user_for_reason()
    data = {
        'reason': reason,
        'tumorSampleBarcode': tumor_sample_barcode
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode('utf-8')

    pubsub_client = pubsub.PublisherClient()
    topic_path = pubsub_client.topic_path(project=project, topic='report.fail')

    message_id = pubsub_client.publish(topic=topic_path, data=encoded_message).result()
    print(f"event emitted! (message id: '{message_id}')")


def _prompt_user_for_reason():
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
    return reason


if __name__ == '__main__':
    main()
