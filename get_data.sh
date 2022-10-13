#!/bin/bash
## Script to download data to recreate results

### Process arguments
while (( "$#" )); do
    case "$1" in
        -h|--help)
            echo "options:"
            echo "-h, --help                show brief help"
            echo "-o, --output-dir=DIR      specify a directory to write downloaded data"
            exit 0
        ;;
        -o)
            shift
            if test $# -gt 0; then
                OUTPUT=$1
            else
                echo "no output dir specified"
                exit 1
            fi
            shift
        ;;
        --output-dir*)
            export OUTPUT=`echo $1 | sed -e 's/^[^=]*=//g'`
            shift
        ;;
        *)
            echo "bad option"
            exit 1
        ;;
    esac
done


##
PARENT="https://bimsbstatic.mdc-berlin.de/akalin/AAkalin_PanelvsMulti/"
wget -P $OUTPUT $PARENT/"data.tar.gz"
tar -xzf $OUTPUT/"data.tar.gz" -C $OUTPUT
rm $OUTPUT/"data.tar.gz"
wget -P $OUTPUT $PARENT/"results.tgz"
tar -xzvf $OUTPUT/"results.tgz" -C $OUTPUT
rm $OUTPUT/"results.tgz"





