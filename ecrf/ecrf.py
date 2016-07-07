from collections import OrderedDict
from lxml import etree
import csv

def field(subject, field_name):
    nodes = subject.xpath('.//ItemData[@ItemOID="%s"]' % field_name)
    return [node.get('Value') for node in nodes]

def code_list(root, field_name):
    code_list_ref = root.xpath('//ItemDef[@OID="%s"]/CodeListRef/@CodeListOID' % field_name)
    if not code_list_ref:
        return {}
    code_list_items = root.xpath('//CodeList[@OID="%s"]/CodeListItem' % code_list_ref[0])
    mapping = {item.get('CodedValue'): item.xpath('.//TranslatedText')[0].text for item in code_list_items}
    return mapping

with open('/tmp/ODM Patientdata.xml') as f:
    root = etree.parse(f).getroot()

    tumor_codes = code_list(root, 'FLD.CARCINOMA.PTUMLOC')

    subjects = root.xpath('//SubjectData')
    subject_data = [OrderedDict([
        ('subject', subject.get('SubjectKey')),
        ('tumor_types', [tumor_codes.get(val) for val in field(subject, 'FLD.CARCINOMA.PTUMLOC')]),
        ('tumor_pct_first_opinions', field(subject, 'FLD.BIOBIOPSIES.FIRSTOP')),
        ('tumor_pct_second_opinions', field(subject, 'FLD.BIOBIOPSIES.SECONDOP')),
    ]) for subject in subjects]

    with open('/tmp/ecrf.csv', 'w') as f:
        writer = csv.DictWriter(f, subject_data[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(subject_data)
