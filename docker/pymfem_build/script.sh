#!/bin/bash
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init -)"
jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --NotebookApp.token='' --NotebookApp.password=''

