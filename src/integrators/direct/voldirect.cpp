/*
 This file is part of Mitsuba, a physically based rendering system.

 Copyright (c) 2007-2014 by Wenzel Jakob and others.

 Mitsuba is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License Version 3
 as published by the Free Software Foundation.

 Mitsuba is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

//#define STATS

#ifdef STATS
static StatsCounter rejectedSamples("Volumetric illumination", "Rejection ratio", EPercentage);
#endif

/*! \plugin{direct}{Direct illumination integrator}
 * \order{1}
 * \parameters{
 *     \parameter{shadingSamples}{\Integer}{This convenience parameter can be
 *         used to set both \code{emitterSamples} and \code{bsdfSamples} at
 *         the same time.
 *     }
 *     \parameter{emitterSamples}{\Integer}{Optional more fine-grained
 *        parameter: specifies the number of samples that should be generated
 *        using the direct illumination strategies implemented by the scene's
 *        emitters\default{set to the value of \code{shadingSamples}}
 *     }
 *     \parameter{bsdfSamples}{\Integer}{Optional more fine-grained
 *        parameter: specifies the number of samples that should be generated
 *        using the BSDF sampling strategies implemented by the scene's
 *        surfaces\default{set to the value of \code{shadingSamples}}
 *     }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See
 *        page~\pageref{sec:strictnormals} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 * }
 * \vspace{-1mm}
 * \renderings{
 *     \medrendering{Only BSDF sampling}{integrator_direct_bsdf}
 *     \medrendering{Only emitter sampling}{integrator_direct_lum}
 *     \medrendering{BSDF and emitter sampling}{integrator_direct_both}
 *     \caption{
 *         \label{fig:integrator-direct}
 *         This plugin implements two different strategies for computing the
 *         direct illumination on surfaces. Both of them are dynamically
 *         combined then obtain a robust rendering algorithm.
 *     }
 * }
 *
 * This integrator implements a direct illumination technique that makes use
 * of \emph{multiple importance sampling}: for each pixel sample, the
 * integrator generates a user-specifiable number of BSDF and emitter
 * samples and combines them using the power heuristic. Usually, the BSDF
 * sampling technique works very well on glossy objects but does badly
 * everywhere else (\subfigref{integrator-direct}{a}), while the opposite
 * is true for the emitter sampling technique
 * (\subfigref{integrator-direct}{b}). By combining these approaches, one
 * can obtain a rendering technique that works well in both cases
 * (\subfigref{integrator-direct}{c}).
 *
 * The number of samples spent on either technique is configurable, hence
 * it is also possible to turn this plugin into an emitter sampling-only
 * or BSDF sampling-only integrator.
 *
 * For best results, combine the direct illumination integrator with the
 * low-discrepancy sample generator (\code{ldsampler}). Generally, the number
 * of pixel samples of the sample generator can be kept relatively
 * low (e.g. \code{sampleCount=4}), whereas the \code{shadingSamples}
 * parameter of this integrator should be increased until the variance in
 * the output renderings is acceptable.
 *
 * \remarks{
 *    \item This integrator does not handle participating media or
 *          indirect illumination.
 * }
 */

class MIVolumetricDirectIntegrator : public SamplingIntegrator
{
public:
	MIVolumetricDirectIntegrator(const Properties &props) :
			SamplingIntegrator(props)
	{
		/* Shorthand notation to set all sample counts */
		size_t shadingSamples = props.getSize("shadingSamples", 1);
		/* Number of samples to take using the emitter sampling technique */
		m_raysamples = props.getSize("raySamples", shadingSamples);
		/* Number of samples to take using the emitter sampling technique */
		m_emitterSamples = props.getSize("emitterSamples", shadingSamples);
		/* Number of samples to take using the BSDF sampling technique */
		m_bsdfSamples = props.getSize("bsdfSamples", shadingSamples);
		/* Number of samples to take using the BSDF sampling technique */
		m_phaseSamples = props.getSize("phaseSamples", shadingSamples);
		/* Maximum number of perfect dielectric of volumetric interfaces */
		m_maxDepth = props.getInteger("maxDepth", -1);
		/* Be strict about potential inconsistencies involving shading normals? */
		m_strictNormals = props.getBoolean("strictNormals", false);
		/* When this flag is set to true, contributions from directly
		 * visible emitters will not be included in the rendered image */
		m_hideEmitters = props.getBoolean("hideEmitters", false);
		/* When this flag is set to false, direct light coming from luminaries
		 * will not be included in the rendered image */
		m_directLight = props.getBoolean("directLight", true);
		/* Choose a random surface point for equiangular sampling, if
		 * false always sample from the luminary center */
		m_sampleEPosition = props.getBoolean("samplePosition", false);

		Assert(m_emitterSamples + m_bsdfSamples > 0);
	}

	/// Unserialize from a binary data stream
	MIVolumetricDirectIntegrator(Stream *stream, InstanceManager *manager) :
			SamplingIntegrator(stream, manager)
	{
		m_raysamples = stream->readSize();
		m_emitterSamples = stream->readSize();
		m_bsdfSamples = stream->readSize();
		m_phaseSamples = stream->readSize();
		m_maxDepth = stream->readInt();
		m_strictNormals = stream->readBool();
		m_hideEmitters = stream->readBool();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const
	{
		SamplingIntegrator::serialize(stream, manager);
		stream->writeSize(m_raysamples);
		stream->writeSize(m_emitterSamples);
		stream->writeSize(m_bsdfSamples);
		stream->writeSize(m_phaseSamples);
		stream->writeInt(m_maxDepth);
		stream->writeBool(m_strictNormals);
		stream->writeBool(m_hideEmitters);
	}

	void configure()
	{
		SamplingIntegrator::configure();

		size_t sumSurf = m_emitterSamples + m_bsdfSamples;
		size_t sumRay = m_emitterSamples + m_phaseSamples;

		m_weightRay = 1.f / (Float) m_raysamples;
		m_weightLum = 1.f / (Float) m_emitterSamples;
		m_weightBSDF = 1.f / (Float) m_bsdfSamples;
		m_weightPhase = 1.f / (Float) m_phaseSamples;

		m_fracBSDF = m_bsdfSamples / (Float) sumSurf;
		m_fracLumSurf = m_emitterSamples / (Float) sumSurf;
		m_fracPhase = m_phaseSamples / (Float) sumRay;
		m_fracLumRay = m_emitterSamples / (Float) sumRay;
	}

	void configureSampler(const Scene *scene, Sampler *sampler)
	{
		SamplingIntegrator::configureSampler(scene, sampler);
		if (m_raysamples > 1)
			sampler->request1DArray(m_raysamples);
		if (m_emitterSamples > 1 || m_raysamples >= 1)
			sampler->request2DArray(m_emitterSamples + m_emitterSamples * m_raysamples);
		if (m_bsdfSamples > 1)
			sampler->request2DArray(m_bsdfSamples);
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const
	{
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		Float eta = 1.0f;

		/* Figure out how many BSDF and direct illumination samples to
		 generate, and where the random numbers should come from */
		bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

		size_t numRaySamples = m_raysamples, numDirectSamples = m_emitterSamples,
				numBSDFSamples = m_bsdfSamples, numPhaseSamples = m_phaseSamples;

		Float fracLumSurf = m_fracLumSurf, fracBSDF = m_fracBSDF, fracLumRay = m_fracLumRay,
				fracPhase = m_fracPhase, weightRay = m_weightRay, weightLum = m_weightLum,
				weightBSDF = m_weightBSDF, weightPhase = m_weightPhase;

		if (rRec.depth > 1 || adaptiveQuery) {
			/* This integrator is used recursively by another integrator.
			 Be less accurate as this sample will not directly be observed. */
			numDirectSamples = 1;
			numRaySamples = std::min(numRaySamples, size_t(1));
			numBSDFSamples = std::min(numBSDFSamples, size_t(1));
			numPhaseSamples = std::min(numPhaseSamples, size_t(1));

			fracBSDF = 0.5f * numBSDFSamples;
			fracLumSurf = 1.f - fracBSDF;
			fracPhase = 0.5f * numPhaseSamples;
			fracLumRay = 1.f - fracPhase;

			weightRay = weightLum = weightBSDF = weightPhase = 1.0f;
		}

		/* Precalculate samples to preserve stratification */
		Float raySample;
		Float* raySampleArray;
		if (numRaySamples > 1) {
			raySampleArray = rRec.sampler->next1DArray(numRaySamples);
		} else {
			raySample = rRec.nextSample1D();
			raySampleArray = &raySample;
		}

		Point2 lumSample, lumRaySample;
		Point2* lumSampleArray;
		if (numDirectSamples > 1 || numRaySamples >= 1) {
			lumSampleArray = rRec.sampler->next2DArray(
					numDirectSamples + numDirectSamples * numRaySamples);
		} else {
			lumSample = rRec.nextSample2D();
			lumSampleArray = &lumSample;
		}
		Point2* lumRaySampleArray = &lumSampleArray[numDirectSamples];

		Point2 bsdfSample;
		Point2* bsdfSampleArray;
		if (numBSDFSamples > 1) {
			bsdfSampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		} else {
			bsdfSample = rRec.nextSample2D();
			bsdfSampleArray = &bsdfSample;
		}

		/* Make sure emitters are invisible */
		RayType oldType = RayType::Default;
		if (m_hideEmitters && (rRec.type & RadianceQueryRecord::ESensorRay)) {
			oldType = ray.type;
			ray.type = ray.type & !RayType::Emmiter;
		}

		/* Perform the first ray intersection (or ignore if the
		 intersection has already been provided). */
		rRec.rayIntersect(ray);

		/* Restore visibility of emitters */
		if (m_hideEmitters && (rRec.type & RadianceQueryRecord::ESensorRay)) {
			ray.type = oldType;
		}

		Spectrum throughput(1.0f);
		bool scattered = false;

		while (rRec.depth < m_maxDepth || m_maxDepth < 0) {
			/* ==================================================================== */
			/*                           Volume sampling                            */
			/* ==================================================================== */
			if (rRec.medium) {
				const PhaseFunction* phase = rRec.medium->getPhaseFunction();

				/* Ray across volume */
				Ray vRay(ray);

				/* Sample volume until nearest intersection */
				if (its.isValid()) {
					vRay.maxt = its.t;
				}

				/* Estimate the single scattering component if this is requested */

				int interactions = m_maxDepth - rRec.depth - 1;
#if EXPONENTIAL
				/* ==================================================================== */
				/*                         Exponential sampling                         */
				/* ==================================================================== */
				for (size_t j = 0; j < numDirectSamples; j++) {
					PositionSamplingRecord pRec;
					scene->sampleEmitterEquiangularNonCenter(vRay, pRec, rRec.medium, interactions,
							lumRaySampleArray[j], raySampleArray[j]);

					const Emitter* emitter = static_cast<const Emitter*>(pRec.object);

					/* Evaluate medium across ray */
					Ray mRay = Ray(vRay, 0.f, distance(pRec.p, vRay.o));

					MediumSamplingRecord mRec;
					rRec.medium->eval(mRay, mRec);
					mRec.t = mRay.maxt;
					mRec.p = pRec.p;

					Spectrum attenuation = mRec.transmittance * mRec.sigmaS;
					DirectSamplingRecord dRec(mRec);

					Spectrum value = emitter->sampleDirect(dRec, lumRaySampleArray[j]);

					value *= scene->evalTransmittance(dRec.ref, false,
							dRec.p, emitter->isOnSurface(), dRec.time, rRec.medium,
							interactions, rRec.sampler);

					if (!value.isZero()) {
						/* Evaluate medium attenuation */
						value *= attenuation;

						/* Evaluate the phase function */
						PhaseFunctionSamplingRecord phRec(mRec, -mRay.d, dRec.d);
						Float phaseVal = phase->eval(phRec);

						if (phaseVal != 0.f) {
							/* Calculate prob. of having sampled that direction using
							 phase function sampling */
							Float phasePdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
							? phase->pdf(phRec) : Float(0.0);

							/* Weight using the power heuristic */
							const Float weight = miWeight(dRec.pdf * fracLum,
									phasePdf * fracBSDF) * weightLum * weightRay;

							Li += throughput * value * phaseVal * weight / pRec.pdf;
						}
					}

				}
#else // EXPONENTIAL
				/* ==================================================================== */
				/*                         Equiangular sampling                         */
				/* ==================================================================== */
				for (size_t i = 0; i < numRaySamples; i++) {
					PositionSamplingRecord pRec;
					if (m_sampleEPosition) {
						scene->sampleEmitterEquiangularNonCenter(vRay, pRec, rRec.medium, interactions,
								rRec.nextSample2D(), raySampleArray[i]);
					} else {
						scene->sampleEmitterEquiangular(vRay, pRec, rRec.medium, interactions,
								raySampleArray[i]);
					}

					if (pRec.pdf == 0)
						continue;

					/* Evaluate medium across ray */
					Ray mRay = Ray(vRay, 0.f, distance(pRec.p, vRay.o));

					MediumSamplingRecord mRec;
					rRec.medium->eval(mRay, mRec);
					mRec.t = mRay.maxt;
					mRec.p = pRec.p;

					Spectrum attenuation = mRec.transmittance * mRec.sigmaS;

					/* ==================================================================== */
					/*                          Luminaire sampling                          */
					/* ==================================================================== */
					DirectSamplingRecord dRec(mRec);

					if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
						const Emitter* emitter = static_cast<const Emitter*>(pRec.object);
#ifdef STATS
						rejectedSamples.incrementBase(numDirectSamples);
#endif
						for (size_t j = 0; j < numDirectSamples; j++) {
							Spectrum value = emitter->sampleDirect(dRec,
									lumRaySampleArray[j + i * numDirectSamples]);

							/* Check rejection sampling */
							size_t rejects = 0;
							while (dRec.measure == EInvalidMeasure && rejects < 100) {
								dRec.measure = ESolidAngle;
								value = emitter->sampleDirect(dRec, rRec.nextSample2D());

								rejects++;
							}
#ifdef STATS
							rejectedSamples += rejects;
#endif
							value *= scene->evalTransmittance(dRec.ref, false, dRec.p,
									emitter->isOnSurface(), dRec.time, rRec.medium,
									interactions, rRec.sampler);

							if (!value.isZero()) {
								/* Evaluate medium attenuation */
								value *= attenuation;

								/* Evaluate the phase function */
								PhaseFunctionSamplingRecord phRec(mRec, -mRay.d, dRec.d);
								Float phaseVal = phase->eval(phRec);

								if (phaseVal != 0.f) {
									/* Calculate prob. of having sampled that direction using
									 phase function sampling */
									Float phasePdf =
											(emitter->isOnSurface()
													&& dRec.measure == ESolidAngle) ?
													phase->pdf(phRec) : Float(0.0);

									/* Weight using the power heuristic */
									const Float weight = miWeight(dRec.pdf * fracLumRay,
											phasePdf * fracPhase) * weightLum * weightRay;

									if (m_directLight || rRec.depth > 1) {
										Li += throughput * value * phaseVal * weight / pRec.pdf;
									}
								}
							}
						}

						/* ==================================================================== */
						/*                         Phase function sampling                      */
						/* ==================================================================== */
						PhaseFunctionSamplingRecord phRec(mRec, -mRay.d);
						dRec = DirectSamplingRecord(mRec);

						for (size_t k = 0; k < numPhaseSamples; k++) {
							Float phasePdf;
							Float phaseVal = phase->sample(phRec, phasePdf, rRec.sampler);

							if (phaseVal == 0.f)
								continue;

							/* Trace a ray in sampled direction */
							Ray phaseRay = Ray(mRec.p, phRec.wo, mRay.time);
							phaseRay.mint = 0.f;

							Spectrum value(0.f);
							Intersection phIts;
							rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
									m_maxDepth - rRec.depth - 1, phaseRay, phIts, dRec, value);

							/* If a luminaire was hit, estimate the local illumination and
							 weight using the power heuristic */
							if (!value.isZero()) {
								/* Evaluate medium attenuation */
								value *= attenuation;

								const Float emitterPdf = scene->pdfEmitterDirect(dRec);

								/* Weight using the power heuristic */;
								const Float weight = miWeight(phasePdf * fracPhase,
										emitterPdf * fracLumRay) * weightPhase * weightRay;
								if (m_directLight || rRec.depth > 1) {
									Li += throughput * value * phaseVal * weight / pRec.pdf;
								}
							}
						}
					}
				}
#endif // EXPONENTIAL
			}

			if (its.isValid()) {
				/* ==================================================================== */
				/*                           Surface sampling                           */
				/* ==================================================================== */
				/* Evaluate medium attenuation (if any) */
				Spectrum surfThroughput(throughput);
				if (rRec.medium) {
					/* Ray across volume */
					Ray vRay(ray);
					vRay.maxt = its.t;

					surfThroughput *= rRec.medium->evalTransmittance(vRay, rRec.sampler);
				}


				/* Possibly include emitted radiance if requested */
				if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
						&& (!m_hideEmitters)) {

					Li += surfThroughput * its.Le(-ray.d);
				}

				const BSDF *bsdf = its.getBSDF(ray);

				/* Estimate the direct illumination if this is requested */
				if ((rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)
						&& (bsdf->getType() & BSDF::ESmooth)) {
					/* ==================================================================== */
					/*                          Luminaire sampling                          */
					/* ==================================================================== */
#ifdef STATS
					rejectedSamples.incrementBase(numDirectSamples);
#endif
					for (size_t i = 0; i < numDirectSamples; i++) {
						int interactions = m_maxDepth - rRec.depth - 1;

						DirectSamplingRecord dRec(its);
						Spectrum value = scene->sampleAttenuatedEmitterDirect(dRec, its,
								rRec.medium, interactions, lumSampleArray[i], rRec.sampler);

						/* Sometimes we sample non-emissive areas and the emitter doesn't register itself */
						if (dRec.object && dRec.measure == EInvalidMeasure) {
							const Emitter* emitter = static_cast<const Emitter*>(dRec.object);

							/* Check rejection sampling */
							size_t rejects = 0;
							while (dRec.measure == EInvalidMeasure && rejects < 100) {
								dRec.measure = ESolidAngle;
								value = emitter->sampleDirect(dRec, rRec.nextSample2D());

								rejects++;
							}
#ifdef STATS
							rejectedSamples += rejects;
#endif
							if (!value.isZero()) {
								dRec.pdf = scene->pdfEmitterDirect(dRec);

								int interactions = m_maxDepth - rRec.depth - 1;
								value *= scene->evalTransmittance(dRec.ref, true, dRec.p,
										emitter->isOnSurface(), dRec.time, rRec.medium,
										interactions, rRec.sampler);
							}
						}

						if (!value.isZero()) {
							const Emitter* emitter = static_cast<const Emitter*>(dRec.object);

							/* Evaluate BSDF * cos(theta) */
							BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
							const Spectrum bsdfVal = bsdf->eval(bRec);

							Float woDotGeoN = dot(its.geoFrame.n, dRec.d);

							/* Prevent light leaks due to the use of shading normals */
							if (!bsdfVal.isZero()
									&& (!m_strictNormals
											|| woDotGeoN * Frame::cosTheta(bRec.wo) > 0.f)) {
								/* Calculate prob. of having generated that direction
								 using BSDF sampling */
								Float bsdfPdf =
										(emitter->isOnSurface() && dRec.measure == ESolidAngle) ?
												bsdf->pdf(bRec) : Float(0.f);

								/* Weight using the power heuristic */
								const Float weight = miWeight(dRec.pdf * fracLumSurf,
										bsdfPdf * fracBSDF) * weightLum;

								if (m_directLight || rRec.depth > 1) {
									Li += surfThroughput * value * bsdfVal * weight;
								}
							}
						}
					}

					/* ==================================================================== */
					/*                            BSDF sampling                             */
					/* ==================================================================== */
					for (size_t j = 0; j < numBSDFSamples; j++) {
						BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);

						Float bsdfPdf;
						Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, bsdfSampleArray[j]);

						if (bsdfVal.isZero())
							continue;

						/* Prevent light leaks due to the use of shading normals */
						const Vector wo = its.toWorld(bRec.wo);
						const Float woDotGeoN = dot(its.geoFrame.n, wo);
						if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
							continue;

						/* Trace a ray in sampled direction */
						DirectSamplingRecord dRec(its);
						Ray bsdfRay = Ray(its.p, wo, ray.time);

						Spectrum value(0.0f);
						Intersection bsdfIts;
						rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
								m_maxDepth - rRec.depth - 1, bsdfRay, bsdfIts, dRec, value);

						/* If a luminaire was hit, estimate the local illumination and
						 weight using the power heuristic */
						if (!value.isZero()) {
							const Float emitterPdf = scene->pdfEmitterDirect(dRec);

							/* Weight using the power heuristic */;
							const Float weight = miWeight(bsdfPdf * fracBSDF,
									emitterPdf * fracLumSurf) * weightBSDF;

							if (m_directLight || rRec.depth > 1) {
								Li += surfThroughput * value * bsdfVal * weight;
							}
						}
					}
				}
			}

			/* ==================================================================== */
			/*                         Multiple scattering                          */
			/* ==================================================================== */
			MediumSamplingRecord mRec;
			if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
				/* ==================================================================== */
				/*                         Phase function sampling                      */
				/* ==================================================================== */
				throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

				const PhaseFunction* phase = mRec.getPhaseFunction();
				DirectSamplingRecord dRec(mRec.p, mRec.time);

				Float phasePdf;
				PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
				Float phaseVal = phase->sample(pRec, phasePdf, rRec.sampler);
				if (phaseVal == 0)
					break;
				throughput *= phaseVal;

				/* Trace a ray in this direction */
				ray = Ray(mRec.p, pRec.wo, ray.time);
				ray.mint = 0.;

				Spectrum value(0.0f);
				rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
						m_maxDepth - rRec.depth - 1, ray, its, dRec, value);

				/* Stop if multiple scattering was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
					break;

				rRec.type = RadianceQueryRecord::ERadianceNoEmission;
			} else {
				/* ==================================================================== */
				/*                            BSDF sampling                             */
				/* ==================================================================== */
				if (rRec.medium)
					throughput *= mRec.transmittance / mRec.pdfFailure;

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return
					 attenuated radiance from a background luminaire */
					if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters)) {
						Spectrum value = throughput * scene->evalEnvironment(ray);
						if (rRec.medium)
							value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
						Li += value;
					}

					break;
				}

				DirectSamplingRecord dRec(mRec.p, mRec.time);

				/* Sample BSDF * cos(theta) */
				BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
				Float bsdfPdf;
				Spectrum bsdfWeight = its.getBSDF(ray)->sample(bRec, bsdfPdf, rRec.nextSample2D());
				if (bsdfWeight.isZero())
					break;

				/* Prevent light leaks due to the use of shading normals */
				const Vector wo = its.toWorld(bRec.wo);
				Float woDotGeoN = dot(its.geoFrame.n, wo);
				if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
					break;

				/* Trace a ray in this direction */
				ray = Ray(its.p, wo, ray.time);

				/* Keep track of the throughput, medium, and relative
				 refractive index along the path */
				throughput *= bsdfWeight;
				eta *= bRec.eta;
				if (its.isMediumTransition())
					rRec.medium = its.getTargetMedium(ray.d);

				/* Handle index-matched medium transitions specially */
				if (bRec.sampledType == BSDF::ENull) {
					if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
						break;
					rRec.type =
							scattered ?
									RadianceQueryRecord::ERadianceNoEmission :
									RadianceQueryRecord::ERadiance;
					scene->rayIntersect(ray, its);
					rRec.depth++;
					continue;
				}

				Spectrum value(0.0f);
				rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
						m_maxDepth - rRec.depth - 1, ray, its, dRec, value);

				/* Stop if indirect illumination was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
					break;

				rRec.type = RadianceQueryRecord::ERadianceNoEmission;
			}

			if (rRec.depth++ >= 4) {
				/* Russian roulette: try to keep path weights equal to one,
				 while accounting for the solid angle compression at refractive
				 index boundaries. Stop with at least some probability to avoid
				 getting stuck (e.g. due to total internal reflection) */
				Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
				if (rRec.nextSample1D() >= q)
					break;
				throughput /= q;
			}

			/* Prepare next loop */
			numDirectSamples = 1;
			numRaySamples = std::min(numRaySamples, size_t(1));
			numBSDFSamples = std::min(numBSDFSamples, size_t(1));
			numPhaseSamples = std::min(numPhaseSamples, size_t(1));

			fracBSDF = 0.5f * numBSDFSamples;
			fracLumSurf = 1.f - fracBSDF;
			fracPhase = 0.5f * numPhaseSamples;
			fracLumRay = 1.f - fracPhase;

			weightRay = weightLum = weightBSDF = weightPhase = 1.0f;

			lumSample = rRec.nextSample2D();
			lumSampleArray = &lumSample;
			if (numRaySamples > 0) {
				raySample = rRec.nextSample1D();
				raySampleArray = &raySample;
				lumRaySample = rRec.nextSample2D();
				lumRaySampleArray = &lumRaySample;
			}
			if (numBSDFSamples > 0) {
				bsdfSample = rRec.nextSample2D();
				bsdfSampleArray = &bsdfSample;
			}

			scattered = true;
		}

		return Li;
	}

	void rayIntersectAndLookForEmitter(const Scene *scene, Sampler *sampler,
		const Medium *medium, int maxInteractions, Ray ray, Intersection &_its,
		DirectSamplingRecord &dRec, Spectrum &value) const
	{
		Intersection its2, *its = &_its;
		Spectrum transmittance(1.0f);
		bool surface = false;
		int interactions = 0;

		while (true) {
			surface = scene->rayIntersect(ray, *its);

			if (medium) {
				transmittance *= medium->evalTransmittance(Ray(ray, 0, its->t), sampler);
			}

			if (surface
					&& (interactions == maxInteractions
							|| !(its->getBSDF()->getType() & BSDF::ENull) || its->isEmitter())) {
				/* Encountered an occluder / light source */
				break;
			}

			if (!surface) {
				break;
			}

			if (transmittance.isZero()) {
				return;
			}

			if (its->isMediumTransition()) {
				medium = its->getTargetMedium(ray.d);
			}

			Vector wo = its->shFrame.toLocal(ray.d);
			BSDFSamplingRecord bRec(*its, -wo, wo, ERadiance);
			bRec.typeMask = BSDF::ENull;
			transmittance *= its->getBSDF()->eval(bRec, EDiscrete);

			ray.o = ray(its->t);
			ray.mint = Epsilon;
			its = &its2;

			if (++interactions > 100) { /// Just a precaution..
				Log(EWarn, "rayIntersectAndLookForEmitter(): round-off error issues?");
				return;
			}
		}

		if (surface) {
			/* Intersected something - check if it was a luminaire */
			if (its->isEmitter()) {
				dRec.setQuery(ray, *its);
				value = transmittance * its->Le(-ray.d);
			}
		} else {
			/* Intersected nothing -- perhaps there is an environment map? */
			const Emitter *env = scene->getEnvironmentEmitter();

			if (env && env->fillDirectSamplingRecord(dRec, ray)) {
				value = transmittance * env->evalEnvironment(RayDifferential(ray));
			}
		}
	}

	inline Float miWeight(Float pdfA, Float pdfB) const
	{
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	std::string toString() const
	{
		std::ostringstream oss;
		oss << "MIVolumetricDirectIntegrator[" << endl << "  emitterSamples = "
				<< m_emitterSamples << "," << endl << "  bsdfSamples = " << m_bsdfSamples << ","
				<< endl << "  maxDepth = " << m_maxDepth << "," << endl << "  strictNormals = "
				<< m_strictNormals << endl << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	size_t m_raysamples;

	size_t m_emitterSamples;
	size_t m_bsdfSamples;
	size_t m_phaseSamples;

	int m_maxDepth;

	Float m_fracBSDF, m_fracPhase, m_fracLumSurf, m_fracLumRay;
	Float m_weightBSDF, m_weightPhase, m_weightLum, m_weightRay;

	bool m_strictNormals;
	bool m_hideEmitters;
	bool m_directLight;
	bool m_sampleEPosition;
};

MTS_IMPLEMENT_CLASS_S(MIVolumetricDirectIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MIVolumetricDirectIntegrator,
		"Direct illumination integrator with volumetric capabilities");
MTS_NAMESPACE_END
