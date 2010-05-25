#!/bin/sh
echo "Start posting T-Coffee latest Vesrion..."
sftp cnotred@tcoffee.org@sftp.tcoffee.org<<EOF
cd mainwebsite_html/Packages/Beta
mput /home/notredame/latest_distributions/T-COFFEE_distribution.tar.gz
EOF

