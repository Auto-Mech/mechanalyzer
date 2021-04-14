# run with: . /path/to/fake-install.sh
export _THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for d in ls $_THIS_DIR/*/; do
    export PATH=$d:$PATH
    export PYTHONPATH=$d:$PYTHONPATH
done

# fake-install external packages as well
. /lcrc/project/CMRP/amech/mechanalyzer/ratefit/external/dsarrfit/debug/fake-install.sh
#. /lcrc/project/CMRP/amech/RateFit/external/troefit/debug/fake-install.sh

