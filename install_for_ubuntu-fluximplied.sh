#!/usr/bin/env bash

printf "Installing dependencies for running fluximplied\n"
sudo apt install -y libxt-dev \
                    libcairo2-dev \
		    libcurl4-openssl-dev \
		    libxml2-dev \
		    libssl-dev
