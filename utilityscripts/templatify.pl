#!/usr/bin/env perl

use 5.16.0;
use strict;
use warnings;

# this script does most of the work for converting the messy, perl-based pipeline jobs into templates
# put all the "print FILEHANDLE" lines into a file and run this on it, redirecting to your template file
# there may be bugs

print <<'HEADER';
#!/usr/bin/env bash
# -*- TT -*-

[% INCLUDE ErrorHandling.tt %]
[% INCLUDE Logging.tt job_name="" main_step= log_name="${step}.log" %]

# you still need to:
# 1) update INCLUDE Logging.tt:
#  a) add job_name (must match template name)
#  b) add "main" step name
#  c) add log file name
#  d) correct success/failure logic
#  e) optionally add sub-steps with start_step "NAME"/finish_step
# 2) template variables inside control structures are potentially double [% %] escaped and should be fixed manually
# 3) inspect manually and convert constructs like \$command, \$mv_command, \$rm_command to that hurts eyes less
# 4) tidy up silly logic, add quoting of paths, formatting into readable multi-line commands etc.

HEADER

while (<>) {
    # beginning/end of line
    s/^\s*print [A-Z_]+ "//;
    s/\\n";\s*$//;

    # remove shebang

    s/^#!.*//;

    # undo escaping

    while (s/^(\s*)\\t/$1    /) {
        next;
    }

    s/\\t/\t/g;
    s/\\n/\n/g;
    s/\\"/"/g;
    s/\\\$/\$/g;
    s/"\.//g;
    s/\."//g;

    # also convert (escaped) backticks to $()

    s/\\`(.+?)(?=\\`)/\$($1)/g;

    # migrate to fail()

    s/echo "ERROR: ([^"]+)"(\s*(>&2|>>.*))?/fail "$1"/;

    # update directory usage

    s#\$opt(?:->)?{OUTPUT_DIR}/(?:\$?[^/\s]+/)*?([^\./\s\$;]+)[\s/"]#[% dirs.$1 %]#g;
    s/\$([a-zA-Z0-9_]+?)_?[dD]ir/[% dirs.$1 %]/g;
    # TODO: fix those underscores
    s#\$opt(?:->)?{OUTPUT_DIR}/#[% dirs.out %]/#g;

    # hash variables first then simple ones

    s/\$([a-zA-Z0-9_]+)(?:->)?{([^}]+)}/[% $1.$2 %]/g;
    s/\$([a-zA-Z0-9_]+)/[% $1 %]/g;

    # conditionals (untested)

    s/^(\s*)if\s*\((.+)\)\s*{\s*$/${1}[% IF $2 -%]/;
    # TODO: FOREACH
    s/^\s*}\s*$/[% END -%]/;

    # paths

    s/opt\.IAP_PATH/opt.OUTPUT_DIR/g;

    say;
}

print <<'FOOTER';

assert_not_empty "[% output_path %]"
success
FOOTER
