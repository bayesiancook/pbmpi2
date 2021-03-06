CC=mpic++
CPPFLAGS= -w -O3 -c
LDFLAGS= -O3
SRCS=  MPIModule.cpp TaxonSet.cpp Tree.cpp Random.cpp SequenceAlignment.cpp CodonSequenceAlignment.cpp \
	NNI.cpp SPR.cpp CCP.cpp BranchSitePath.cpp \
	MultiGeneMPIModule.cpp MultiGeneRateProcess.cpp MultiGeneBranchProcess.cpp MultiGeneProfileProcess.cpp MultiGeneSBDPProfileProcess.cpp \
	MultiGenePhyloProcess.cpp MultiGeneRASCATGTRSBDPGammaPhyloProcess.cpp MultiGeneRASCATSBDPGammaPhyloProcess.cpp MultiGeneAASubSelRASCATSBDPGammaPhyloProcess.cpp \
	StateSpace.cpp CodonStateSpace.cpp ZippedSequenceAlignment.cpp SubMatrix.cpp \
	ZipSubMatrix.cpp GTRSubMatrix.cpp CodonSubMatrix.cpp linalg.cpp Chrono.cpp BranchProcess.cpp \
	GammaBranchProcess.cpp RateProcess.cpp DGamRateProcess.cpp ProfileProcess.cpp DirichletProfileProcess.cpp DirichletMixtureProfileProcess.cpp \
	SiteSpecificProfileProcess.cpp MatrixSiteSpecificProfileProcess.cpp GeneralPathSuffStatMatrixSiteSpecificProfileProcess.cpp \
	CovMatrix.cpp MVNProfileProcess.cpp MVNSiteSpecificProfileProcess.cpp AACodonMutSelMVNSiteSpecificPhyloProcess.cpp \
	GTRProfileProcess.cpp GTRDPProfileProcess.cpp GTRSBDPProfileProcess.cpp GTRFiniteProfileProcess.cpp ExpoConjugateGTRProfileProcess.cpp \
	MixtureProfileProcess.cpp ZipMatrixProfileProcess.cpp \
	PoissonMixtureProfileProcess.cpp MatrixMixtureProfileProcess.cpp MultipleMatrixMixtureProfileProcess.cpp \
	GTRMixtureProfileProcess.cpp ExpoConjugateGTRMixtureProfileProcess.cpp \
	GeneralPathSuffStatMatrixMixtureProfileProcess.cpp GeneralPathSuffStatMultipleMatrixMixtureProfileProcess.cpp \
	DPProfileProcess.cpp SBDPProfileProcess.cpp \
	PoissonDPProfileProcess.cpp PoissonSBDPProfileProcess.cpp \
	FiniteProfileProcess.cpp PoissonFiniteProfileProcess.cpp \
	CodonMutSelProfileProcess.cpp \
	AACodonMutSelSBDPOmegaSBDPPhyloProcess.cpp AACodonMutSelSBDPOmegaFinitePhyloProcess.cpp SBDPOmegaProcess.cpp AACodonMutSelSiteOmegaFinitePhyloProcess.cpp \
	AACodonMutSelSiteOmegaSBDPPhyloProcess.cpp SiteOmegaProcess.cpp AACodonMutSelSiteSBDPPhyloProcess.cpp AACodonMutSelFiniteOmegaSBDPPhyloProcess.cpp \
	OmegaProcess.cpp MixtureOmegaProcess.cpp SingleOmegaProcess.cpp FiniteOmegaProcess.cpp AACodonMutSelProfileProcess.cpp \
	SubstitutionProcess.cpp Propagate.cpp Parsimony.cpp \
	PoissonSubstitutionProcess.cpp MatrixSubstitutionProcess.cpp \
	ExpoConjugateGTRSubstitutionProcess.cpp GeneralPathSuffStatMatrixSubstitutionProcess.cpp ExpoConjugateGTRPhyloProcess.cpp \
	ZipGeneralPathSuffStatMatrixMixtureProfileProcess.cpp ZipMatrixPhyloProcess.cpp ZipExpoConjugateGTRPhyloProcess.cpp ReadPhyloProcess.cpp \
	ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess.cpp ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess.cpp \
	SMC.cpp PhyloProcess.cpp PoissonPhyloProcess.cpp MatrixPhyloProcess.cpp \
	GeneralPathSuffStatMatrixPhyloProcess.cpp \
	RASCATGammaPhyloProcess.cpp RASCATSBDPGammaPhyloProcess.cpp \
	RASCATFiniteGammaPhyloProcess.cpp RASCATGTRFiniteGammaPhyloProcess.cpp \
	RASCATGTRDPGammaPhyloProcess.cpp RASCATGTRSBDPGammaPhyloProcess.cpp \
	GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess.cpp \
	GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess.cpp \
	AACodonMutSelFinitePhyloProcess.cpp AACodonMutSelFiniteOmegaFinitePhyloProcess.cpp \
	CodonMutSelFinitePhyloProcess.cpp\
	CodonMutSelSBDPPhyloProcess.cpp \
	AACodonMutSelSiteMatrixMixtureProfileProcess.cpp AACodonMutSelSBDPPhyloProcess.cpp \
	AASubSelSubMatrix.cpp AASubSelMixtureProfileProcess.cpp AASubSelRASCATSBDPGammaPhyloProcess.cpp \
	PartitionDGamRateProcess.cpp PartitionGTRProfileProcess.cpp PartitionGTRMixtureProfileProcess.cpp PartitionGTRFiniteProfileProcess.cpp \
	ExpoConjugatePartitionGTRProfileProcess.cpp ExpoConjugatePartitionGTRMixtureProfileProcess.cpp \
	ExpoConjugatePartitionGTRSubstitutionProcess.cpp ExpoConjugatePartitionGTRPhyloProcess.cpp \
	RASPARTGTRGammaPhyloProcess.cpp \
	GammaRateProcess.cpp PoissonSiteSpecificProfileProcess.cpp RASIIDDirichletGammaPhyloProcess.cpp IIDDirichletIIDGammaPhyloProcess.cpp \
	Bipartition.cpp BipartitionList.cpp Consensus.cpp TaxaParameters.cpp PBTree.cpp TreeList.cpp PolyNode.cpp correl.cpp correlation.cpp \
	BP2Stat.cpp BP2util.cpp


OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= pb_mpi readpb_mpi logo tracecomp simupartcodonrev simupartaamix 
#ALL= logo satplot pb_mpi readpb_mpi geneempfreq reroot hi stickbreaking drawhindex gettreelength simupartaamix simupartcodonrev partitionjackknife simucodonrev diversity mccpmcmc ccpmcmc bpcomp bpcomp2 tracecomp
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

$(PROGSDIR)/logo: MakeSequenceLogo.o $(OBJS)
	$(CC) MakeSequenceLogo.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/stickbreaking: StickBreaking.o $(OBJS)
	$(CC) StickBreaking.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/drawhindex: DrawHindex.o $(OBJS)
	$(CC) DrawHindex.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/satplot: SaturationPlot.o $(OBJS)
	$(CC) SaturationPlot.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/hi: HI.o $(OBJS)
	$(CC) HI.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/geneempfreq: GeneEmpiricalFreq.o $(OBJS)
	$(CC) GeneEmpiricalFreq.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/reroot: Reroot.o $(OBJS)
	$(CC) Reroot.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/gettreelength : GetTreeLength.o $(OBJS)
	$(CC) GetTreeLength.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/simupartaamix: SimuPartitionAAMix.o $(OBJS)
	$(CC) SimuPartitionAAMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/simupartcodonrev: SimuPartitionCodonRev.o $(OBJS)
	$(CC) SimuPartitionCodonRev.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/partitionjackknife: PartitionJackknife.o $(OBJS)
	$(CC) PartitionJackknife.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/mccpmcmc: MultiCCPMCMC.o $(OBJS)
	$(CC) MultiCCPMCMC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/ccpmcmc: CCPMCMC.o $(OBJS)
	$(CC) CCPMCMC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/diversity: Diversity.o $(OBJS)
	$(CC) Diversity.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/missfrac: MissingFractionPerGene.o $(OBJS)
	$(CC) MissingFractionPerGene.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/jackknife: JackKnife.o $(OBJS)
	$(CC) JackKnife.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/simucodonrev: SimuCodonRev.o $(OBJS)
	$(CC) SimuCodonRev.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/simucodon: SimuCodon.o $(OBJS)
	$(CC) SimuCodon.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/readpb_mpi: ReadPB.o $(OBJS)
	$(CC) ReadPB.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/bpcomp2: BPcomp2.o $(OBJS)
	$(CC) BPcomp2.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/pb_mpi: PB.o $(OBJS)
	$(CC) PB.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/tracecomp: TraceComp.o $(OBJS)
	$(CC) TraceComp.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/bpcomp: BPCompare.o $(OBJS)
	$(CC) BPCompare.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/woconst: woconst.o $(OBJS)
	$(CC) woconst.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

