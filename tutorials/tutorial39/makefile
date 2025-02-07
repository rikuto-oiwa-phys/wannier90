outfiles := CRASH \
bands.out \
bandsx.out \
cwf_parameters.dat \
nscf.out \
out/ \
pw2wan.out \
p_w_vs_e.dat \
scf.out \
silicon.amn \
silicon.band \
silicon_band.dat \
silicon_band.pdf \
silicon.band.gnu \
silicon_band.gnu \
silicon_band.kpt \
silicon_band.labelinfo.dat \
silicon_u* \
silicon.chk \
silicon.eig \
silicon.mmn \
silicon.nnkp \
silicon.werr \
silicon.wout \
UNK* \
*xsf

infiles := silicon.bands \
silicon.bandsx \
silicon.nscf \
silicon.pw2wan \
silicon.scf \
silicon.win \
silicon_bandsdiff.gnu \

default:
	@echo "Usage: make clean"

clean:
	@echo "Cleaning up..."
	rm -r $(outfiles)

# Usually, it is mpirun. If not set, run in serial.
# MPI_CMD = mpirun -n 8

test:
	@echo "Testing all calculations"
	$(MPI_CMD) pw.x < silicon.scf > scf.out
	$(MPI_CMD) pw.x < silicon.bands > bands.out
	$(MPI_CMD) bands.x < silicon.bandsx > bandsx.out
	$(MPI_CMD) pw.x < silicon.nscf > nscf.out
	wannier90.x -pp silicon
	$(MPI_CMD) pw2wannier90.x < silicon.pw2wan > pw2wan.out
	../../utility/fit_cwf_parameters.py silicon
	MUMAX=`grep 'cwf_mu_max' cwf_parameters.dat | awk '{print $$3}'` ; \
	MUMIN=`grep 'cwf_mu_min' cwf_parameters.dat | awk '{print $$3}'` ; \
	SIGMAX=`grep 'cwf_sigma_max' cwf_parameters.dat | awk '{print $$3}'` ; \
	SIGMIN=`grep 'cwf_sigma_min' cwf_parameters.dat | awk '{print $$3}'` ; \
    sed "/cwf_mu_max/s/^cwf_mu_max.*$/cwf_mu_max     = ${MUMAX}/g" silicon.win > temp.win && mv temp.win silicon.win ; \
	sed "/cwf_mu_min/s/^cwf_mu_min.*$/cwf_mu_min     = ${MUMIN}/g" silicon.win > temp.win && mv temp.win silicon.win ; \
	sed "/cwf_sigma_max/s/^cwf_sigma_max.*$/cwf_sigma_max  = ${SIGMAX}/g" silicon.win > temp.win && mv temp.win silicon.win ; \
	sed "/cwf_sigma_min/s/^cwf_sigma_min.*$/cwf_sigma_min  = ${SIGMIN}/g" silicon.win > temp.win && mv temp.win silicon.win ; \
        $(MPI_CMD) wannier90.x silicon
	./silicon_bandsdiff.gnu
	@echo "All calculations finished successfully"
