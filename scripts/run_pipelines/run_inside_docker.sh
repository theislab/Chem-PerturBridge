SCRIPT_DIR="$(dirname "$0")"

#Create environment
if [ ! -d ~/.venv ]; then
    echo "create .venv"
    chmod u+x $SCRIPT_DIR/set_env.sh
    $SCRIPT_DIR/set_env.sh
fi

#./pipelines/sciplex/sciplex_pseudobulking_subsample.sh
#./pipelines/sciplex/sciplex_pseudobulking.sh
./pipelines/tahoe/tahoe_pseudobulking_subsample.sh
#./pipelines/tahoe/tahoe_pseudobulking.sh
