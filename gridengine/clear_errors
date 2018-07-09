SGE_ROOT=/opt/sge; export SGE_ROOT
for i in $(/opt/sge/bin/linux-x64/qstat|grep Eqw|cut -f 1 -d " ");do /opt/sge/bin/linux-x64/qmod -c $i;done
