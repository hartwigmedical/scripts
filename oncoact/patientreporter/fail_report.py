import json
import argparse
from google.cloud import pubsub
from cli_util import perform_prod_test


def main():
    """
    This script is used to easily emit a fail report command/event to the reporting pipeline.
    """
    parser = argparse.ArgumentParser(description="Emit a QC fail event to the reporting pipeline.")
    parser.add_argument('tumor_sample_barcode')
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    perform_prod_test(args.profile)

    tumor_sample_barcode = args.tumor_sample_barcode
    project = 'hmf-pipeline-development'  # TODO make dynamic
    assemble_and_emit_event(tumor_sample_barcode, project)


def assemble_and_emit_event(tumor_sample_barcode, project):
    """
    Gather all the input information needed to emit the report failed event.

    This method will prompt the user for the relevant information.

    :param tumor_sample_barcode: the tumor sample barcode of the report to emit the event for
    :param project: the gcp project to emit this event to.
    """
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

    topic_path = client.topic_path(project, 'report.fail')
    res = emit_fail_report(client, topic_path, reason, tumor_sample_barcode)

    print(f"event emitted! (message id: '{res}')")
    exit(0)


def emit_fail_report(client: pubsub.PublisherClient,
                     topic_path: str,
                     tumor_sample_barcode: str,
                     reason: str) -> str:
    """
    Emits a report failed event.

    :param client: the pubsub client.
    :param topic_path: the topic path (which includes both project and topic).
    :param tumor_sample_barcode: the tumor sample barcode to emit the event for.
    :param reason: the reason for failure.
    :return: the id of the emitted event
    """
    data = {
        'reason': reason,
        'tumorSampleBarcode': tumor_sample_barcode
    }
    json_data = json.dumps(data)
    encoded_message = json_data.encode('utf-8')
    return client.publish(topic=topic_path, data=encoded_message).result()


if __name__ == '__main__':
    main()
