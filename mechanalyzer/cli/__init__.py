import click
from mechanalyzer.cli import ste_mech


@click.group()
def main():
    """MechAnalyzer CLI"""
    pass


@main.command()
def expand():
    """Expand stereochemistry for a mechanism
    """
    ste_mech.main()


@main.command()
def greetme():
    """Hello world function, for CLI testing purposes"""
    print("Hello, world!")
