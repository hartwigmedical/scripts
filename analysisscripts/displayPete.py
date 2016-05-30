__author__ = 'dhowe'

from enum import Enum#, unique
from IPython.core.display import HTML, display
import numpy as np


#@unique
class FormatType(Enum):
    billions_2 = 1
    billions_1 = 5
    thousands_1 = 2
    ones_0 = 3
    ones_1 = 15
    ones_2 = 4
    ones_3 = 12
    bp_2 = 7
    pc_2 = 8
    pc_0 = 9
    nf = 6
    bool = 10
    millions_2 = 11
    millions_1 = 13
    pc_1 = 14
    ones_0_nc = 16


def formatter_pete(x, format_type):
    green = '#006600'
    red = '#CC0000'
    black = '#000000'
    no_comma = False
    try:
        if np.isnan(x):
            return ""
    except TypeError:
        pass
    if format_type == FormatType.nf:
        return x
    elif format_type == FormatType.bool:
        return '<font color="{}"><div align = right>{}'.format(green if x else red, x)
    elif format_type == FormatType.billions_2:
        modified_x = x/1000000000
        num_decimals = 2
        suffix = 'bn'
    elif format_type == FormatType.billions_1:
        modified_x = x/1000000000
        num_decimals = 1
        suffix = 'bn'
    elif format_type == FormatType.ones_0:
        modified_x = x
        num_decimals = 0
        suffix = ''
    elif format_type == FormatType.ones_1:
        modified_x = x
        num_decimals = 1
        suffix = ''
    elif format_type == FormatType.ones_2:
        modified_x = x
        num_decimals = 2
        suffix = ''
    elif format_type == FormatType.ones_3:
        modified_x = x
        num_decimals = 3
        suffix = ''
    elif format_type == FormatType.thousands_1:
        modified_x = x/1000
        num_decimals = 1
        suffix = 'k'
    elif format_type == FormatType.bp_2:
        modified_x = x
        num_decimals = 2
        suffix = 'bp'
    elif format_type == FormatType.pc_2:
        modified_x = x*100
        num_decimals = 2
        suffix = '%'
    elif format_type == FormatType.pc_0:
        modified_x = x*100
        num_decimals = 0
        suffix = '%'
    elif format_type == FormatType.pc_1:
        modified_x = x*100
        num_decimals = 1
        suffix = '%'
    elif format_type == FormatType.millions_1:
        modified_x = x/1000000
        num_decimals = 1
        suffix = 'mn'
    elif format_type == FormatType.millions_2:
        modified_x = x/1000000
        num_decimals = 2
        suffix = 'mn'
    elif format_type == FormatType.ones_0_nc:
        modified_x = x
        num_decimals = 0
        suffix = ''
        no_comma = True
    color = red if modified_x <= -5*(10**(-num_decimals-1)) else black
    if no_comma:
        return '<font color="{{}}"><div align = right>{{:.{}f}}{{}}'.format(num_decimals).format(color, modified_x,
                                                                                                  suffix)
    else:
        return '<font color="{{}}"><div align = right>{{:,.{}f}}{{}}'.format(num_decimals).format(color, modified_x,
                                                                                                  suffix)

formatter_dict = dict()
#for format_type in FormatType:
#    formatter_dict[format_type] = lambda x, format_type=format_type: formatter_pete(x, format_type)


def display_for_pete(df, format_type_list, return_html_string=False):
    # myFormatter = lambda x: '<font color="#FF0000"><div align = right>{:,.0f}'.format(x) if x<0
    # else '<font color="#000000"><div align = right>{:,.0f}'.format(x)
    # ## This one formats with commas and red for negative numbers

    my_format_dict = {"nothing": formatter_pete}

    for i, c in enumerate(df.columns):
        my_format_dict[c] = formatter_dict[format_type_list[i]]
    if return_html_string:
        return df.to_html(formatters=my_format_dict, escape=False)
    else:
        display(HTML(df.to_html(formatters=my_format_dict, escape=False)))

# Not exactly sure what the below does but it may be important so adding it in.
style = """
<style>
.container { width:90% !important; }
</style>
"""
HTML(style)


def display_for_pete_format_rows(df, format_type_list):
    rows = df.index
    values = df.values
    headings = df.keys()
    formatted_values = [[formatter_pete(x, format_type_list[i]) for x in row] for i, row in enumerate(values)]

    html = """<table border="1" class="dataframe"> \n
                  <thead> \n
                    <tr style="text-align: right;"> \n
                      <th></th> \n"""
    for heading in headings:
        html += '<th>' + str(heading) + '</th> \n'
    html += """</tr> \n
            </thead> \n
            <tbody> \n
                <tr> \n"""
    for row, value_list in zip(rows, formatted_values):
        html += '<tr> \n'
        html += '<th>' + str(row) + '</th> \n'
        for value in value_list:
            html += '<td>' + str(value) + '</td> \n'
        html += '</tr> \n'
    html += """ </tbody> \n
                </table>"""
    display(HTML(html))