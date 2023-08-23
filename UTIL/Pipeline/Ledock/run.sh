#!/bin/sh

cp /opt/scripts/ledock.py ./
python ledock.py --config config.in --run_type ${1}
