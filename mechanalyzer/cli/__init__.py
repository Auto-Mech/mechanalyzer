import click

from mechanalyzer.cli import run_sort, ste_mech


@click.group()
def main():
    """MechAnalyzer CLI"""
    pass


@main.command()
@click.option(
    "-m",
    "--mech",
    default="mechanism.dat",
    show_default=True,
    help="Input mechanism file name",
)
@click.option(
    "-s",
    "--spc",
    default="species.csv",
    show_default=True,
    help="Input species file name",
)
@click.option(
    "-t",
    "--therm",
    default="therm.dat",
    show_default=True,
    help="Input thermo file name",
)
@click.option(
    "-i",
    "--sort",
    default="sort.dat",
    show_default=True,
    help="Input sort file name",
)
@click.option(
    "-o",
    "--outmech",
    default="outmech.dat",
    show_default=True,
    help="Output mechanism file name",
)
@click.option(
    "-c",
    "--outspc",
    default="outspc.csv",
    show_default=True,
    help="Output species file name",
)
@click.option(
    "-g",
    "--outgroups",
    default="pes_groups.dat",
    show_default=True,
    help="Output PES groups file name",
)
def sort(
    mech: str = "mechanism.dat",
    spc: str = "species.csv",
    therm: str = "therm.dat",
    sort: str = "sort.dat",
    outmech: str = "outmech.dat",
    outspc: str = "outspc.csv",
    outgroups: str = "pes_groups.dat",
):
    """Sort the reactions in a mechanism"""
    run_sort.main(
        mech=mech,
        spc=spc,
        therm=therm,
        sort=sort,
        outmech=outmech,
        outspc=outspc,
        outgroups=outgroups,
    )


@main.command()
def expand():
    """Expand stereochemistry for a mechanism"""
    ste_mech.main()
