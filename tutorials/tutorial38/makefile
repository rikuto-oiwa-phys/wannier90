outfiles := CRASH \
bands.out \
bandsx.out \
nscf.out \
out/ \
pw2wan.out \
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
silicon_bandsdiff.gnu

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
	$(MPI_CMD) wannier90.x silicon
	./silicon_bandsdiff.gnu
	@echo "All calculations finished successfully"
