#include "integrator.h"

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO
    for (int i = 0; i < particleSystem.getCubeCount(); i += 1) {
        for (int j = 0; j < particleSystem.cubes[i].getParticleNum(); j += 1) {
            Eigen::Vector3f force = particleSystem.cubes[i].getParticle(j).getForce();
            Eigen::Vector3f velocity = particleSystem.cubes[i].getParticle(j).getVelocity();
            Eigen::Vector3f position = particleSystem.cubes[i].getParticle(j).getPosition();
            float mass = particleSystem.cubes[i].getParticle(j).getMass();
            float deltatime = particleSystem.deltaTime;

            // newpos = origin_pos + vt
            position += velocity * deltatime;
            // f = ma, v = at, v = ft/m
            velocity += force * deltatime / mass;
            // set new pos and v
            particleSystem.cubes[i].getParticle(j).setVelocity(velocity);
            particleSystem.cubes[i].getParticle(j).setPosition(position);

            /*printf("%.4lf ", particleSystem.cubes[i].getParticle(0).getForce().dot(Eigen::Vector3f(1, 0, 0)));
            printf("%.4lf ", particleSystem.cubes[i].getParticle(0).getForce().dot(Eigen::Vector3f(0, 1, 0)));
            printf("%.4lf ", particleSystem.cubes[i].getParticle(0).getForce().dot(Eigen::Vector3f(0, 0, 1)));
            printf("\n");*/
        }
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO
    for (int i = 0; i < particleSystem.getCubeCount(); i += 1)
    {
        particleSystem.computeCubeForce(particleSystem.cubes[i]);
        for (int j = 0; j < particleSystem.cubes[i].getParticleNum(); j += 1)
        {
            Eigen::Vector3f force = particleSystem.cubes[i].getParticle(j).getForce();
            Eigen::Vector3f velocity = particleSystem.cubes[i].getParticle(j).getVelocity();
            Eigen::Vector3f position = particleSystem.cubes[i].getParticle(j).getPosition();
            float mass = particleSystem.cubes[i].getParticle(j).getMass();
            float delta_time = particleSystem.deltaTime;

            //  f = ma, v = v0 + at, v = v0 + f / m * t
            velocity = force / mass * delta_time;
            //  newpos = oripos + vt
            position += velocity * delta_time;

            particleSystem.cubes[i].getParticle(j).setVelocity(velocity);
            particleSystem.cubes[i].getParticle(j).setPosition(position);
        }
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO
    // For midpoint euler, the deltaTime passed in is correct.
    // But this deltaTime is for a full step.
    // So you may need to adjust it before computing, but don't forget to restore original value.
    for (int i = 0; i < particleSystem.getCubeCount(); i += 1)
    {
        for (int j = 0; j < particleSystem.cubes[i].getParticleNum(); j += 1)
        {
            Eigen::Vector3f force = particleSystem.cubes[i].getParticle(j).getForce();
            Eigen::Vector3f velocity = particleSystem.cubes[i].getParticle(j).getVelocity();
            Eigen::Vector3f position = particleSystem.cubes[i].getParticle(j).getPosition();
            float mass = particleSystem.cubes[i].getParticle(j).getMass();
            float delta_time = particleSystem.deltaTime;

            //velocity += delta_time / 2 * particleSystem.cubes[i].getParticle(j).getAcceleration();
            velocity += delta_time / 2 * force / mass;
            position += delta_time * velocity;
            //  f = ma, v = v + at, v = v + f / m * t
            velocity = particleSystem.cubes[i].getParticle(j).getVelocity() + delta_time * force / mass ;

            particleSystem.cubes[i].getParticle(j).setVelocity(velocity);
            particleSystem.cubes[i].getParticle(j).setPosition(position);
        }
    }
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    // TODO
    // StateStep struct is just a hint, you can use whatever you want.
    for (int i = 0; i < particleSystem.getCubeCount(); i += 1)
    {
        for (int j = 0; j < particleSystem.cubes[i].getParticleNum(); j += 1)
        {
            Eigen::Vector3f force = particleSystem.cubes[i].getParticle(j).getForce();
            Eigen::Vector3f velocity = particleSystem.cubes[i].getParticle(j).getVelocity();
            Eigen::Vector3f position = particleSystem.cubes[i].getParticle(j).getPosition();
            float mass = particleSystem.cubes[i].getParticle(j).getMass();
            float delta_time = particleSystem.deltaTime;

            StateStep s;
            s.deltaVel = Eigen::Vector3f(0, 0, 0);
            s.deltaPos = Eigen::Vector3f(0, 0, 0);
            Eigen::Vector3f V = Eigen::Vector3f(0, 0, 0);
            Eigen::Vector3f P = Eigen::Vector3f(0, 0, 0);
            // k1
            s.deltaVel += velocity + force / mass * delta_time;
            s.deltaPos = position + s.deltaVel * delta_time;
            V += s.deltaVel;
            P += s.deltaPos;
            // k2
            s.deltaVel += s.deltaVel * delta_time/2;
            s.deltaPos = position + delta_time * s.deltaVel;
            V += s.deltaVel * 2;
            P += s.deltaPos * 2;
            // k3
            s.deltaVel += s.deltaVel * delta_time / 2;
            s.deltaPos = position + delta_time * s.deltaVel;
            V += s.deltaVel * 2;
            P += s.deltaPos * 2;
            // k4
            s.deltaVel += s.deltaVel * delta_time;
            s.deltaPos = position + delta_time * s.deltaVel;
            V += s.deltaVel;
            P += s.deltaPos;

            particleSystem.cubes[i].getParticle(j).setVelocity(V/6);
            particleSystem.cubes[i].getParticle(j).setPosition(P/6);
        }
    }
}
}  // namespace simulation
