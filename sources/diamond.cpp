
class ProfileProcess	{

};

class MixtureProfileProcess	{

	MixtureProfileProcess() {}
	virtual ~MixtureProfileProcess() {}

	virtual double ProfileSuffStatLogProb(int k) = 0;
};

class MatrixMixtureProfileProcess : public virtual MixtureProfileProcess	{

};

