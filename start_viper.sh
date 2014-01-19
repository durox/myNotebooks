#!/bin/sh

ipython1x notebook --pylab inline --no-browser --port 8888 &
ipython=$!

ipython1x viper.py

kill $ipython
