#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 1.2
  // Using BSDF::reflect(), implement sample_f for a mirror surface
  reflect(wo, wi);
  *pdf = 1;
  return reflectance / abs_cos_theta(*wi);

  //return Spectrum();
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  float distribution;
  float exponent = (float) -pow((sin_theta(h)/cos_theta(h)), 2) / (alpha*alpha);
  distribution = (float) (exp(exponent) / (PI * (alpha*alpha) * pow(cos_theta(h),4)));
  return distribution;


  //return std::pow(cos_theta(h), 100.0);;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.

  Spectrum eta_and_k = (eta * eta) + (k * k);
  float cos_squared = (float) (cos_theta(wi) * cos_theta(wi));
  Spectrum two_eta_cos = 2*eta*cos_theta(wi);
  Spectrum rs = (eta_and_k - two_eta_cos + cos_squared) / (eta_and_k + two_eta_cos + cos_squared);
  Spectrum rp = (eta_and_k * cos_squared - two_eta_cos + 1) / (eta_and_k * cos_squared + two_eta_cos + 1);
  return (rs+rp) / 2;

  //return Spectrum();
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here

  Vector3D n = Vector3D(0,0,1);
  float wi_check = (float) dot(n,wi);
  float wo_check = (float) dot(n,wo);
  if (wi_check > 0 && wo_check > 0) {
    Vector3D h = wo + wi;
    h.normalize();
    return (F(wi) * G(wo, wi) * D(h)) / (4 * (dot(n, wo)) * dot(n, wi));
  }
  return Spectrum();
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

  float theta, phi;
  Vector2D samples = sampler.get_sample();
  float sample1 = (float)samples[0];
  float sample2 = (float)samples[1];
  theta = atan(sqrt(-(alpha*alpha)*log(1-sample1)));
  phi = (float)(2*PI*sample2);
  Vector3D h = Vector3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  *wi = -wo + 2 * dot(wo, h) * h;
  if (wi->z < 0) {
    *pdf = 0;
    return Spectrum();
  }
  float exponent = (float) (-pow(tan(theta), 2) / (alpha*alpha));
  float pdf_theta = (float) ((2*sin(theta) / ((alpha*alpha)*pow(cos(theta),3))) * exp(exponent));
  float pdf_phi = (float) (1 / (2 * PI));
  float pdf_sampling_h = pdf_theta*pdf_phi/sin(theta);
  *pdf = (float) (pdf_sampling_h/(4*dot(*wi,h)));
  return MicrofacetBSDF::f(wo, *wi);

//  *wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
//  return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.
  if (!refract(wo, wi, ior)) {
    reflect(wo, wi);
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
  } else {
    float r0, r;
    r0 = ((1-ior)/(1+ior)) * ((1-ior)/(1+ior));
    r = (float) (r0 + (1-r0) * pow((1-abs_cos_theta(wo)),5));
    if (coin_flip(r)) {
      reflect(wo, wi);
      *pdf = r;
      return r * reflectance / abs_cos_theta(*wi);
    } else {
      refract(wo, wi, ior);
      *pdf = 1-r;
      float n;
      if (wo.z < 0) {
        n = ior;
      } else {
        n = 1/ior;
      }
      return (1-r) * transmittance / abs_cos_theta(*wi) / (n*n);
    }
  }

}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0], -wo[1], wo[2]);

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  float n;
  float condition;

  if (wo.z < 0) {
    n = ior;
    condition = (float)(1-(n*n)*(1-(wo.z*wo.z)));
    if (condition < 0) {
      return false;
    }
    *wi = Vector3D(-n*wo.x, -n*wo.y, sqrt(condition));
    return true;
  } else {
    n = 1/ior;
    condition = (float)(1-(n*n)*(1-(wo.z*wo.z)));
    if (condition < 0) {
      return false;
    }
    *wi = Vector3D(-n*wo.x, -n*wo.y, -sqrt(condition));
    return true;
  }
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
