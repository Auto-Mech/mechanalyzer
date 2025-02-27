import click
import numpy as np

from mechanalyzer.cli import (sort,
                              prompt,
                              pssa,
                              ste_mech, 
                              compare_rates, 
                              #compare_thermo
                              )


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
    "--sortopts",
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
def sortmech(
    mech: str = "mechanism.dat",
    spc: str = "species.csv",
    therm: str = "therm.dat",
    sortopts: str = "sort.dat",
    outmech: str = "outmech.dat",
    outspc: str = "outspc.csv",
    outgroups: str = "pes_groups.dat",
    ):
    """Sort the reactions in a mechanism"""
    sort.main(
        mech=mech,
        spc=spc,
        therm=therm,
        sortopts=sortopts,
        outmech=outmech,
        outspc=outspc,
        outgroups=outgroups,
    )

@main.command()
@click.option(
    "-m",
    "--mechs_yaml",
    default="mechs.yaml",
    show_default=True,
    help="yaml file, set up like:\\mech1:\\\trate_file: 'rate.ckin'\\\ttherm_file: 'therm.ckin'\\\tspecies_csv: 'species.csv'\\mech2: ...",
)
@click.option(
    "-o",
    "--plot_fname",
    default="rate_plot.pdf",
    show_default=True,
    help="name of output pdf of plots"
)
@click.option(
    "-f",
    "--out_txt_fname",
    default="ordering.txt",
    show_default=True,
    help="name of output text filename"
)
@click.option(
    "-d",
    "--job_path",
    default="",
    show_default=True,
    help="directory for input/output files"
)
@click.option(
    "-t",
    "--temps_lst",
    default=None,
    type=lambda s: [float(temp) for temp in s.split(',')],
    show_default=True,
    help="array of temperatures, None defaults to 500--1500"
)
@click.option(
    "-p",
    "--pressures",
    default=None,
    type=lambda s: [float(press) for press in s.split(',')],
    show_default=True,
    help="array of pressures to plot, None defaults to (1, 10, 100)"
)
@click.option(
    "-s",
    "--sort_method",
    default="ratios",
    show_default=True,
    help="Sort the pdf plots by the ratio of the differences with 'ratios', or not at all with None"
)
@click.option(
    "-r",
    "--rev_rates",
    default=True,
    show_default=True,
    help="reverse rates of remaining mechanisms to make them match the direction in the first"
)
@click.option(
    "-l",
    "--remove_loners",
    default=True,
    show_default=True,
    help="True only plots rates that are in ALL mechanisms, False plots ALL rates in all mechanism"
)
def compare_mechanisms(
    mechs_yaml: str = 'mechs.yaml',
    plot_fname: str = 'rate_plot.pdf',
    out_txt_fname: str = 'ordering.txt',
    job_path: str = '',
    temps_lst: list = [],
    pressures: list = [],
    sort_method: str = None,
    rev_rates: bool = True,
    remove_loners: bool = True,
):
    """Compare the rate constants in a mechanism"""
    compare_rates.main(
        mechs_yaml=mechs_yaml,
        plot_fname=plot_fname,
        out_txt_fname=out_txt_fname,
        job_path=job_path,
        temps_lst=temps_lst,
        pressures=pressures,
        sort_method=sort_method,
        rev_rates=rev_rates,
        remove_loners=remove_loners
    )

@main.command()
@click.option(
    "-m",
    "--mechs_yaml",
    default="mechs.yaml",
    show_default=True,
    help="yaml file, set up like:\\mech1:\\\ttherm_file: 'therm.ckin'\\\tspecies_csv: 'species.csv'\\mech2: ...",
)
@click.option(
    "-o",
    "--plot_fname",
    default="thermo_plot.pdf",
    show_default=True,
    help="name of output pdf of plots"
)
@click.option(
    "-f",
    "--out_txt_fname",
    default="ordering.txt",
    show_default=True,
    help="name of output text filename"
)
@click.option(
    "-d",
    "--job_path",
    default=".",
    show_default=True,
    help="directory for input/output files"
)
@click.option(
    "-t",
    "--temps_lst",
    default=None,
    type=lambda s: [float(temp) for temp in s.split(',')],
    show_default=True,
    help="array of temperatures, None defaults to 500--1500"
)
@click.option(
    "-s",
    "--sort_method",
    default="lnq",
    show_default=True,
    help="Sort the pdf plots by the max difference in enthalpy ('h'), entropy ('s'), Gibbs ('g'), c_p ('cp'),\\natural log of the partition function ('lnq'), or not at all (None)"
)
@click.option(
    "-st",
    "--sort_temp",
    default=None,
    show_default=True,
    help="If sorting, specifies the temp at which to sort. None (default) specifies to sort by maximum difference."
)
@click.option(
    "-l",
    "--remove_loners",
    default=True,
    show_default=True,
    help="True only plots species that are in ALL mechanisms, False plots ALL species in all mechanism"
)
@click.option(
    "-p",
    "--print_missing",
    default=True,
    show_default=True,
    help="True prints a warning for any species that are not in the species.csv file"
)
def compare_thermo(
    mechs_yaml: str = 'mechs.yaml',
    plot_fname: str = 'thermo_plot.pdf',
    out_txt_fname: str = 'ordering.txt',
    job_path: str = '.',
    temps_lst: list = [],
    sort_method: str = None,
    sort_temp: float = None,
    rev_rates: bool = True,
    remove_loners: bool = True,
    print_missing: bool = True
):
    """Compare the thermo properties in a mechanism"""
    compare_thermo.main(
        mechs_yaml=mechs_yaml,
        plot_fname=plot_fname,
        out_txt_fname=out_txt_fname,
        job_path=job_path,
        temps_lst=temps_lst,
        sort_method=sort_method,
        sort_temp=sort_temp,
        remove_loners=remove_loners,
        print_missing=print_missing
    )

@main.command()
def expand():
    """Expand stereochemistry for a mechanism"""
    ste_mech.main()


@main.command()
@click.option(
    "-f",
    "--flds",
    default= ['',],
    show_default=True,
    required=True,
    multiple=True,
    help="directories with mess files; pass each directory as a single -f argument",
)
@click.option(
    "-i",
    "--messinput",
    default="mess.inp",
    show_default=True,
    help="mess input",
)
@click.option(
    "-o",
    "--messoutput",
    default="mess.out",
    show_default=True,
    help="mess rate output",
)
@click.option(
    "-ke",
    "--outputmicro",
    default="ke.out",
    show_default=True,
    help="mess microcanonical rate output",
)
@click.option(
    "-l",
    "--log",
    default="mess.log",
    show_default=True,
    help="mess log file",
)
@click.option(
    "-m",
    "--model",
    default="rovib_dos",
    show_default=True,
    help="Model to compute prompt branching fractions",
)
@click.option(
    "-b",
    "--bfthresh",
    default=0.1,
    show_default=True,
    help="keep reactions if contributing more than this threshold",
)
@click.option(
    "-fit",
    "--fitmethod",
    default="plog",
    show_default=True,
    help="Fitting method",
)
@click.option(
    "-or",
    "--outputrates",
    default="rates_prompt.txt",
    show_default=True,
    help="Output prompt rates file name",
)
def promptcalc(
    flds: list = ['',],
    messinput: str = 'mess.inp',
    messoutput: str = 'mess.out',
    outputmicro: str = 'ke.out',
    log: str = 'mess.log',
    model: str = 'rovib_dos',
    bfthresh: float = 0.1,
    fitmethod: str = 'plog',
    outputrates: str = 'rates_prompt.txt'
):
    """Compute prompt effects from mess output files"""
    prompt.main(
        flds=flds,
        messinput=messinput,
        messoutput=messoutput,
        outputmicro=outputmicro,
        log=log,
        model=model,
        bfthresh=bfthresh,
        fitmethod=fitmethod,
        outputrates=outputrates
    )

@main.command()
@click.option(
    "-i",
    "--startmech",
    default='kin.CKI',
    show_default=True,
    required=True,
    help="input mechanism in chemkin format",
)
@click.option(
    "-s",
    "--pssa_spcs",
    default=['',],
    show_default=True,
    required=True,
    multiple=True,
    help="list of species to apply pssa to; provide each -s argument separately",
)
@click.option(
    "-T",
    "--temprange",
    default=[500, 2000],
    show_default=True,
    multiple=True,
    help="temperature range for rate constant calculation; provide 2 limits as -T arguments separately",
)
@click.option(
    "-P",
    "--pvect",
    default=[1.0,],
    show_default=True,
    multiple=True,
    type=float,
    help="pressure vector for rate constant calculation; provide pressure values as multiple -P arguments separately",
)
@click.option(
    "-b",
    "--bfthresh",
    default=1e-4,
    show_default=True,
    help="keep reactions if contributing more than this threshold",
)
@click.option(
    "-t",
    "--thermofile",
    default=None,
    show_default=True,
    help="thermochemistry file to compute backward rate constants",
)
@click.option(
    "-or",
    "--outputrates",
    default="rates_prompt.txt",
    show_default=True,
    help="Output prompt rates file name",
)
def runpssa(
    startmech: str='kin.CKI',
    pssa_spcs: list=['',],
    temprange: list = [500, 2000],
    pvect: list =[1.0,],
    bfthresh: float = 1e-4,
    thermofile: str = None,
    outputrates: str = 'pssa_rates.txt',
):
    """Sort the reactions in a mechanism"""
    pssa.main(
    startmech=startmech,
    pssa_spcs=pssa_spcs,
    temprange=temprange,
    pvect=pvect,
    bfthresh=bfthresh,
    thermofile=thermofile,
    outputrates=outputrates,
    )