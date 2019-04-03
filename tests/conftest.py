import os

import click
import click.testing
import pytest
from typing import TextIO, List


class _FileTester(object):
    def __init__(self, cli: click.Command = None, args: List[str] = None):
        self.runner = click.testing.CliRunner(mix_stderr=False)
        self.cli = cli
        self.result = self.invoke(cli, args) \
            if cli is not None and args is not None else None

    @classmethod
    def output_lines(cls, file, sort: bool = False):
        lines = [l.strip() for l in file.readlines()]
        return sorted(lines) if sort else lines

    @classmethod
    def output_string(cls, file, sort: bool = False):
        return '\n'.join(str(line) for line in cls.output_lines(file, sort))

    def invoke(self, cli: click.Command, args: List[str],
               catch_exceptions=False, verbose=True) -> click.testing.Result:
        """Run run a command with the given arguments. Essentially
        a convenience function for click.testing.CliRunner.invoke.
        Exceptions will be raised unless catch_exceptions is set to True.

        :param cli: a click.command function to be run
        :param args: a list of command line arguments to invoke it with
        :param catch_exceptions: if set to True (not recommended) the
            exceptions will be caught by click and limited exception
            information will be emitted
        :param verbose: show the command line
        :return: the result of running the command
        :rtype click.testing.Result:
        """
        if verbose:
            click.echo(f"\nRunning: {cli.name} {' '.join(args)}", err=True)
        result = self.runner.invoke(cli, args,
                                    catch_exceptions=catch_exceptions)

        return result

    def assert_correct_output(self,
                              correct_output_filename: str,
                              out_file: TextIO = None,
                              verbose: bool = False,
                              test_order: bool = True) -> click.testing.Result:
        """Test the output of the command against the correct output file.
        The command is run only if no command
        Specify a file name if output should of the command should be evaluated
        against a different file than the standard output. First the command
        is asserted to have zero exit status and then to have the correct
        output.

        :param bool test_order: if set to False, ignore the order of the output
            lines. Compare the sorted lines.
        :param str correct_output_filename: a filename to compare against
        :param TextIO out_file: Output file to be compared against.
            This should only be specified when the output to be compared is not
            standard output.
        :param bool verbose: write stderr to test.
        :return: The result object from running the command.
        :rtype: click.testing.Result
        """
        if self.result is None:
            raise ValueError("Command was not run yet")

        if verbose:
            if len(self.result.stderr):
                click.echo("\n" + self.result.stderr)

        output = self.output_lines(open(correct_output_filename),
                                   sort=not test_order)

        if out_file is None:
            lines = self.result.output.strip().split('\n')
            correct_output = lines if test_order else sorted(lines)
            assert correct_output == output
        else:
            assert self.output_lines(out_file, sort=not test_order) == output
        return self.result

    @classmethod
    def _result_stderr(cls, result: click.testing.Result):
        """Workaround: click.testing.Result fails to properly check stderr
        separation in its property and throws exceptions unnecessarily."""
        if len(result.stderr_bytes):
            return result.stderr
        else:
            return u""

    def assert_exit_code(self,
                         expected_exit_code: int = 0,
                         verbose=True) -> click.testing.Result:
        """Test the exit code against the expected exit code when running a
            command.

        :param expected_exit_code: expected exit status of the code
        :return:
        """
        if self.result is None:
            raise ValueError("Command was not run yet")

        assert self.result.exit_code == expected_exit_code, \
            f"Exit code for {self.cli.name} was {self.result.exit_code}, " \
            f"not {expected_exit_code}.\nCommand: {self.cli.name} " \
            f"{' '.join(self.cli.params)}\noutput:\n" \
            f"{self._result_stderr(self.result)}\n"

        return self.result


@pytest.fixture
def FileTester():
    return _FileTester


@pytest.fixture
def output_dir():
    return os.path.join(os.path.dirname(__file__), 'output')


@pytest.fixture
def fasta_dir():
    return os.path.join(os.path.dirname(__file__), 'fasta')
