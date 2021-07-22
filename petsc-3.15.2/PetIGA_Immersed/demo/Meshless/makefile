



TARGETS = \
FSIzhang \
FSIstick2 \
FSIdamage \
Fluid5 \
Vacuum3D

ALL: ${TARGETS}
clean::
	-@${RM} ${TARGETS}

CFLAGS    = -w
FFLAGS    =
CPPFLAGS  = -std=c++11
FPPFLAGS  =
LOCDIR    = demo/Meshless/
EXAMPLESC =
EXAMPLESF =
MANSEC    = IGA

topdir := $(shell cd ../.. && pwd)
PETIGA_DIR ?= $(topdir)
include ${PETIGA_DIR}/conf/petigavariables
include ${PETIGA_DIR}/conf/petigarules

Column: FSI3Dcolumn.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<
Phase: PhaseBlast2.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<
Box: test.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB} -lboost_graph_parallel -lboost_mpi -lboost_system -lboost_serialization
	${RM} -f $<
Fluid: Fluidnew.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<
Fluid2: RDX_immersed.o chkopts
		${CLINKER} -o $@ $< ${PETIGA_LIB}
		${RM} -f $<
VacuumNSNI_SinglePlate: VacuumSC_NSNI_SinglePlate.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<
VacuumnoNSNI_SinglePlate: VacuumSC_noNSNI_SinglePlate.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<
VacuumNSNI_2plate: VacuumSC_NSNI_2Plate.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<
VacuumnoNSNI_2plate: VacuumSC_noNSNI_2Plate.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<
Bar_Vibration_NSNI: Bar_NSNI.o chkopts
	${CLINKER} -o $@ $< ${PETIGA_LIB}
	${RM} -f $<

include ${PETIGA_DIR}/conf/petigatest
