#!/usr/bin/env python3

import sys
import os
import argparse
import shutil
import subprocess
from .__version__ import __version__
import textwrap as _textwrap
from .align import align


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
DIM = '\033[2m'
SUPPRESS = '==SUPPRESS=='
OPTIONAL = '?'
ZERO_OR_MORE = '*'
ONE_OR_MORE = '+'
PARSER = 'A...'
REMAINDER = '...'
_UNRECOGNIZED_ARGS_ATTR = '_unrecognized_args'

def main():
    args = parse_args(sys.argv[1:])
    align(args)


def parse_args(args):
    description = BOLD+'gapmm2: gapped alignment with minimap2.'+END_FORMATTING+\
        ' Performs minimap2/mappy alignment with splice options and refines terminal alignments with edlib. Output is PAF format.'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)
    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('reference', help='reference genome (FASTA)')
    required_args.add_argument('query', help='transcipts in FASTA or FASTQ')

    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument('-o', '--out', help='output in PAF format (default: stdout)', metavar='')
    optional_args.add_argument('-t', '--threads', type=int, default=3, help='number of threads to use with minimap2 (default: 3)', metavar='')
    optional_args.add_argument('-m', '--min-mapq', dest='min_mapq', default=1, help='minimum map quality value', metavar='')
    optional_args.add_argument('-i', '--max-intron', dest='max_intron', type=int, default=500, help='max intron length, controls terminal search space', metavar='')
    optional_args.add_argument('-d', '--debug', action='store_true', help='write some debug info to stderr')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='GFFtk v' + __version__,
                           help="Show program's version number and exit")

    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)

## formatting clases for argparse below
class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if a command is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """
    def error(self, message):
        if len(sys.argv) == 2:  # if a command was given but nothing else
            self.print_help(file=sys.stderr)
            sys.exit(2)
        else:
            super().error(message)

class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        self.colours = get_colours_from_tput()
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if 'default' not in help_text.lower():
                help_text += ' (default: {})'.format(action.default)
            elif 'default: DEFAULT' in help_text:
                help_text = help_text.replace('default: DEFAULT',
                                              'default: {}'.format(action.default))
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = BOLD + heading + END_FORMATTING
        super().start_section(heading)

    def _split_lines(self, text, width):
        """
        Override this method to add special behaviour for help texts that start with:
          'R|' - loop text one option per line
        """
        if text.startswith('R|'):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            for line in text_lines:
                if len(line) <= width:
                    wrapped_text_lines.append(line)
                else:
                    wrap_column = 2
                    line_parts = line.split(', ')
                    join = ','
                    current_line = line_parts[0]
                    for part in line_parts[1:]:
                        if len(current_line) + len(join) + 1 + len(part) <= width:
                            current_line += join + ' ' + part
                        else:
                            wrapped_text_lines.append(current_line + join)
                            current_line = ' ' * wrap_column + part
                    wrapped_text_lines.append(current_line)
            return wrapped_text_lines
        elif text.startswith('C|'):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            for line in text_lines:
                if len(line) <= width:
                    wrapped_text_lines.append(line)
                else:
                    wrap_column = 2
                    line_parts = line.split(', ')
                    join = ','
                    current_line = line_parts[0]
                    for part in line_parts[1:]:
                        if len(current_line) + len(join) + 1 + len(part) <= width:
                            current_line += join + ' ' + part
                        else:
                            wrapped_text_lines.append(current_line + join)
                            current_line = ' ' * wrap_column + part
                    wrapped_text_lines.append(current_line)
            return wrapped_text_lines
        else:
            return argparse.HelpFormatter._split_lines(self, text, width)

    def _fill_text(self, text, width, indent):
        if text.startswith('R|'):
            return ''.join(indent + line for line in text[2:].splitlines(keepends=True))
        elif text.startswith('C|'):
            return ''.join(line for line in text[2:].splitlines(keepends=True))
        else:
            text = self._whitespace_matcher.sub(' ', text).strip()
            paragraphs = text.split('|n')
            multiline_text = ''
            for paragraph in paragraphs:
                formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n'
                multiline_text = multiline_text + formatted_paragraph
            return multiline_text

    def _format_args(self, action, default_metavar):
        get_metavar = self._metavar_formatter(action, default_metavar)
        if action.nargs is None:
            result = '%s' % get_metavar(1)
        elif action.nargs == OPTIONAL:
            result = '[%s]' % get_metavar(1)
        elif action.nargs == ZERO_OR_MORE:
            metavar = get_metavar(1)
            if len(metavar) == 2:
                result = '[%s [%s ...]]' % metavar
            else:
                result = '[%s ...]' % metavar
        elif action.nargs == ONE_OR_MORE:
            result = '%s' % get_metavar(1)
        elif action.nargs == REMAINDER:
            result = '...'
        elif action.nargs == PARSER:
            result = '%s ...' % get_metavar(1)
        elif action.nargs == SUPPRESS:
            result = ''
        else:
            try:
                formats = ['%s' for _ in range(action.nargs)]
            except TypeError:
                raise ValueError("invalid nargs value") from None
            result = ' '.join(formats) % get_metavar(action.nargs)
        return result


def get_colours_from_tput():
    try:
        return int(subprocess.check_output(['tput', 'colors']).decode().strip())
    except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
        return 1


if __name__ == '__main__':
    main()