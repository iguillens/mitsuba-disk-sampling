#include <iostream>

#include "DiskSampler.hpp"
#include "vector3.hpp"

int main ()
{
	DiskSampler<double> sampler = 
		DiskSampler<double>(SphericalEllipseSampler<double>::Type::POLARTAB);

	v3l::vector3<double> o = v3l::vector3<double>(0.0, 0.0, 0.0);
	v3l::vector3<double> cd = v3l::vector3<double>(0.0, 3.0, 3.0);
	v3l::vector3<double> nd = v3l::vector3<double>(0.0, -1.0, 0.0);
	double r = 3.0;

	DiskSamplingRecord<double> samplingRecord = sampler.createRecord(o, cd, nd, r);

	double u = 0.15, v = 0.5;

	double pdf = 0.0;
	v3l::vector3<double> p = samplingRecord.sample(u, v, pdf);

	std::cout << "p=" << p << std::endl;
	std::cout << "pdf=" << pdf << std::endl;

	std::cout << "****************" << std::endl;

	pdf = samplingRecord.samplePdf(p);
	std::cout << "pdf=" << pdf << std::endl;
}
