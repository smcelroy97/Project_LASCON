// #define DEBUG 1
/*
 *  adex_gamma_E_ml.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Generated from NESTML at time: 2024-01-24 19:08:12.837180
**/

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "adex_gamma_E_ml.h"

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
nest::RecordablesMap<adex_gamma_E_ml> adex_gamma_E_ml::recordablesMap_;
namespace nest
{

  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
template <> void RecordablesMap<adex_gamma_E_ml>::create()
  {
    // add state variables to recordables map
   insert_(adex_gamma_E_ml_names::_V_m, &adex_gamma_E_ml::get_V_m);
   insert_(adex_gamma_E_ml_names::_w, &adex_gamma_E_ml::get_w);
   insert_(adex_gamma_E_ml_names::_g_AMPA__X__AMPA, &adex_gamma_E_ml::get_g_AMPA__X__AMPA);
   insert_(adex_gamma_E_ml_names::_g_NMDA__X__NMDA, &adex_gamma_E_ml::get_g_NMDA__X__NMDA);
   insert_(adex_gamma_E_ml_names::_g_NMDA__DOLLAR__X__NMDA, &adex_gamma_E_ml::get_g_NMDA__DOLLAR__X__NMDA);
   insert_(adex_gamma_E_ml_names::_g_GABA__X__GABA, &adex_gamma_E_ml::get_g_GABA__X__GABA);
    // add recordable inline expressions to recordables map
	insert_(adex_gamma_E_ml_names::_I_AMPA, &adex_gamma_E_ml::get_I_AMPA);
	insert_(adex_gamma_E_ml_names::_I_GABA, &adex_gamma_E_ml::get_I_GABA);
	insert_(adex_gamma_E_ml_names::_I_NMDA, &adex_gamma_E_ml::get_I_NMDA);

    // Add vector variables  
  }
}
std::vector< std::tuple< int, int > > adex_gamma_E_ml::rport_to_nestml_buffer_idx =
{
  { adex_gamma_E_ml::AMPA, adex_gamma_E_ml::PORT_NOT_AVAILABLE },
  { adex_gamma_E_ml::NMDA, adex_gamma_E_ml::PORT_NOT_AVAILABLE },
  { adex_gamma_E_ml::GABA, adex_gamma_E_ml::PORT_NOT_AVAILABLE },
};

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of adex_gamma_E_ml's constructor.
// ---------------------------------------------------------------------------

adex_gamma_E_ml::Parameters_::Parameters_()
{
}

adex_gamma_E_ml::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

adex_gamma_E_ml::Buffers_::Buffers_(adex_gamma_E_ml &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

adex_gamma_E_ml::Buffers_::Buffers_(const Buffers_ &, adex_gamma_E_ml &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

adex_gamma_E_ml::adex_gamma_E_ml():ArchivingNode(), P_(), S_(), B_(*this)
{
  init_state_internal_();
  recordablesMap_.create();
  pre_run_hook();
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

adex_gamma_E_ml::adex_gamma_E_ml(const adex_gamma_E_ml& __n):
  ArchivingNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this) {

  // copy parameter struct P_
  P_.C_m = __n.P_.C_m;
  P_.t_ref = __n.P_.t_ref;
  P_.V_reset = __n.P_.V_reset;
  P_.g_L = __n.P_.g_L;
  P_.E_L = __n.P_.E_L;
  P_.a = __n.P_.a;
  P_.b = __n.P_.b;
  P_.Delta_T = __n.P_.Delta_T;
  P_.tau_w = __n.P_.tau_w;
  P_.V_th = __n.P_.V_th;
  P_.V_peak = __n.P_.V_peak;
  P_.E_e = __n.P_.E_e;
  P_.E_i = __n.P_.E_i;
  P_.tau_decay_AMPA = __n.P_.tau_decay_AMPA;
  P_.tau_decay_GABA = __n.P_.tau_decay_GABA;
  P_.Mg = __n.P_.Mg;
  P_.tau_syn_rise_NMDA = __n.P_.tau_syn_rise_NMDA;
  P_.tau_syn_decay_NMDA = __n.P_.tau_syn_decay_NMDA;
  P_.I_e = __n.P_.I_e;

  // copy state struct S_
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::w] = __n.S_.ode_state[State_::w];
  S_.r = __n.S_.r;
  S_.ode_state[State_::g_AMPA__X__AMPA] = __n.S_.ode_state[State_::g_AMPA__X__AMPA];
  S_.ode_state[State_::g_NMDA__X__NMDA] = __n.S_.ode_state[State_::g_NMDA__X__NMDA];
  S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA] = __n.S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA];
  S_.ode_state[State_::g_GABA__X__GABA] = __n.S_.ode_state[State_::g_GABA__X__GABA];

  // copy internals V_
  V_.t_peak_NMDA = __n.V_.t_peak_NMDA;
  V_.__h = __n.V_.__h;
  V_.g_NMDA_const = __n.V_.g_NMDA_const;
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  V_.__P__g_AMPA__X__AMPA__g_AMPA__X__AMPA = __n.V_.__P__g_AMPA__X__AMPA__g_AMPA__X__AMPA;
  V_.__P__g_NMDA__X__NMDA__g_NMDA__X__NMDA = __n.V_.__P__g_NMDA__X__NMDA__g_NMDA__X__NMDA;
  V_.__P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA = __n.V_.__P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA;
  V_.__P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA = __n.V_.__P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA;
  V_.__P__g_GABA__X__GABA__g_GABA__X__GABA = __n.V_.__P__g_GABA__X__GABA__g_GABA__X__GABA;
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

adex_gamma_E_ml::~adex_gamma_E_ml()
{
  // GSL structs may not have been allocated, so we need to protect destruction

  if (B_.__s)
  {
    gsl_odeiv_step_free( B_.__s );
  }

  if (B_.__c)
  {
    gsl_odeiv_control_free( B_.__c );
  }

  if (B_.__e)
  {
    gsl_odeiv_evolve_free( B_.__e );
  }
}

// ---------------------------------------------------------------------------
//   Node initialization functions
// ---------------------------------------------------------------------------
void adex_gamma_E_ml::calibrate_time( const nest::TimeConverter& tc )
{
  LOG( nest::M_WARNING,
    "adex_gamma_E_ml",
    "Simulation resolution has changed. Internal state and parameters of the model have been reset!" );

  init_state_internal_();
}
void adex_gamma_E_ml::init_state_internal_()
{
#ifdef DEBUG
  std::cout << "adex_gamma_E_ml::init_state_internal_()" << std::endl;
#endif

  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  // use a default "good enough" value for the absolute error. It can be adjusted via `node.set()`
  P_.__gsl_error_tol = 1e-3;
  // initial values for parameters
    

    P_.C_m = 150.0; // as pF
    

    P_.t_ref = 5.0; // as ms
    

    P_.V_reset = (-65.0); // as mV
    

    P_.g_L = 10.0; // as nS
    

    P_.E_L = (-65.0); // as mV
    

    P_.a = 4; // as nS
    

    P_.b = 20; // as pA
    

    P_.Delta_T = 2.0; // as mV
    

    P_.tau_w = 500; // as ms
    

    P_.V_th = (-40); // as mV
    

    P_.V_peak = 30; // as mV
    

    P_.E_e = 0; // as mV
    

    P_.E_i = (-80); // as mV
    

    P_.tau_decay_AMPA = 1.5; // as ms
    

    P_.tau_decay_GABA = 7.5; // as ms
    

    P_.Mg = 1; // as real
    

    P_.tau_syn_rise_NMDA = 2; // as ms
    

    P_.tau_syn_decay_NMDA = 200; // as ms
    

    P_.I_e = 0; // as pA

  recompute_internal_variables();
  // initial values for state variables
    

    S_.ode_state[State_::V_m] = (-65); // as mV
    

    S_.ode_state[State_::w] = 0; // as pA
    

    S_.r = 0; // as integer
    

    S_.ode_state[State_::g_AMPA__X__AMPA] = 0; // as real
    

    S_.ode_state[State_::g_NMDA__X__NMDA] = 0; // as real
    

    S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA] = 0; // as real
    

    S_.ode_state[State_::g_GABA__X__GABA] = 0; // as real
}

void adex_gamma_E_ml::init_buffers_()
{
#ifdef DEBUG
  std::cout << "adex_gamma_E_ml::init_buffers_()" << std::endl;
#endif
  // spike input buffers
  get_spike_inputs_().clear();
  get_spike_inputs_grid_sum_().clear();

  // continuous time input buffers  

  get_I_stim().clear();
  B_.I_stim_grid_sum_ = 0;

  B_.logger_.reset();



  if ( not B_.__s )
  {
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.__s );
  }

  if ( not B_.__c )
  {
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( not B_.__e )
  {
    B_.__e = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = adex_gamma_E_ml_dynamics;
  B_.__sys.jacobian = nullptr;
  B_.__sys.dimension = State_::STATE_VEC_SIZE;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void adex_gamma_E_ml::recompute_internal_variables(bool exclude_timestep) {
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  if (exclude_timestep) {    
      

      V_.t_peak_NMDA = P_.tau_syn_decay_NMDA * P_.tau_syn_rise_NMDA * std::log(P_.tau_syn_decay_NMDA / P_.tau_syn_rise_NMDA) / (P_.tau_syn_decay_NMDA - P_.tau_syn_rise_NMDA); // as real
      

      V_.g_NMDA_const = 1 / (std::exp((-V_.t_peak_NMDA) / P_.tau_syn_decay_NMDA) - std::exp((-V_.t_peak_NMDA) / P_.tau_syn_rise_NMDA)); // as real
      

      V_.RefractoryCounts = nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps(); // as integer
      

      V_.__P__g_AMPA__X__AMPA__g_AMPA__X__AMPA = std::exp((-V_.__h) / P_.tau_decay_AMPA); // as real
      

      V_.__P__g_NMDA__X__NMDA__g_NMDA__X__NMDA = std::exp((-V_.__h) / P_.tau_syn_rise_NMDA); // as real
      

      V_.__P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA = (-P_.tau_syn_decay_NMDA) * P_.tau_syn_rise_NMDA * std::exp((-V_.__h) / P_.tau_syn_rise_NMDA) / (P_.tau_syn_decay_NMDA - P_.tau_syn_rise_NMDA) + P_.tau_syn_decay_NMDA * P_.tau_syn_rise_NMDA * std::exp((-V_.__h) / P_.tau_syn_decay_NMDA) / (P_.tau_syn_decay_NMDA - P_.tau_syn_rise_NMDA); // as real
      

      V_.__P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA = std::exp((-V_.__h) / P_.tau_syn_decay_NMDA); // as real
      

      V_.__P__g_GABA__X__GABA__g_GABA__X__GABA = std::exp((-V_.__h) / P_.tau_decay_GABA); // as real
  }
  else {    
      

      V_.t_peak_NMDA = P_.tau_syn_decay_NMDA * P_.tau_syn_rise_NMDA * std::log(P_.tau_syn_decay_NMDA / P_.tau_syn_rise_NMDA) / (P_.tau_syn_decay_NMDA - P_.tau_syn_rise_NMDA); // as real
      

      V_.__h = __resolution; // as ms
      

      V_.g_NMDA_const = 1 / (std::exp((-V_.t_peak_NMDA) / P_.tau_syn_decay_NMDA) - std::exp((-V_.t_peak_NMDA) / P_.tau_syn_rise_NMDA)); // as real
      

      V_.RefractoryCounts = nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps(); // as integer
      

      V_.__P__g_AMPA__X__AMPA__g_AMPA__X__AMPA = std::exp((-V_.__h) / P_.tau_decay_AMPA); // as real
      

      V_.__P__g_NMDA__X__NMDA__g_NMDA__X__NMDA = std::exp((-V_.__h) / P_.tau_syn_rise_NMDA); // as real
      

      V_.__P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA = (-P_.tau_syn_decay_NMDA) * P_.tau_syn_rise_NMDA * std::exp((-V_.__h) / P_.tau_syn_rise_NMDA) / (P_.tau_syn_decay_NMDA - P_.tau_syn_rise_NMDA) + P_.tau_syn_decay_NMDA * P_.tau_syn_rise_NMDA * std::exp((-V_.__h) / P_.tau_syn_decay_NMDA) / (P_.tau_syn_decay_NMDA - P_.tau_syn_rise_NMDA); // as real
      

      V_.__P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA = std::exp((-V_.__h) / P_.tau_syn_decay_NMDA); // as real
      

      V_.__P__g_GABA__X__GABA__g_GABA__X__GABA = std::exp((-V_.__h) / P_.tau_decay_GABA); // as real
  }
}
void adex_gamma_E_ml::pre_run_hook() {
  B_.logger_.init();

  // parameters might have changed -- recompute internals
  recompute_internal_variables();

  // buffers B_
  B_.spike_inputs_.resize(NUM_SPIKE_RECEPTORS);
  B_.spike_inputs_grid_sum_.resize(NUM_SPIKE_RECEPTORS);
}

// ---------------------------------------------------------------------------
//   Update and spike handling functions
// ---------------------------------------------------------------------------

extern "C" inline int adex_gamma_E_ml_dynamics(double __time, const double ode_state[], double f[], void* pnode)
{
  typedef adex_gamma_E_ml::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const adex_gamma_E_ml& node = *( reinterpret_cast< adex_gamma_E_ml* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  f[State_::V_m] = node.P_.E_e * (1.0 * ode_state[State_::g_NMDA__X__NMDA] / (0.280112044817927 * node.P_.C_m * node.P_.Mg * std::exp((-0.062) * std::min(ode_state[State_::V_m], node.P_.V_peak)) + node.P_.C_m) + 1.0 * ode_state[State_::g_AMPA__X__AMPA] / node.P_.C_m) - 1.0 * ode_state[State_::g_NMDA__X__NMDA] * std::min(ode_state[State_::V_m], node.P_.V_peak) / (0.280112044817927 * node.P_.C_m * node.P_.Mg * std::exp((-0.062) * std::min(ode_state[State_::V_m], node.P_.V_peak)) + node.P_.C_m) + node.P_.E_L * node.P_.g_L / node.P_.C_m + 1.0 * node.P_.E_i * ode_state[State_::g_GABA__X__GABA] / node.P_.C_m + node.P_.I_e / node.P_.C_m + (node.P_.Delta_T * node.P_.g_L * std::exp(((-node.P_.V_th) + std::min(ode_state[State_::V_m], node.P_.V_peak)) / node.P_.Delta_T) - 1.0 * ode_state[State_::g_AMPA__X__AMPA] * std::min(ode_state[State_::V_m], node.P_.V_peak) - 1.0 * ode_state[State_::g_GABA__X__GABA] * std::min(ode_state[State_::V_m], node.P_.V_peak) - node.P_.g_L * std::min(ode_state[State_::V_m], node.P_.V_peak) - ode_state[State_::w]) / node.P_.C_m;
  f[State_::w] = node.P_.a * ((-node.P_.E_L) / node.P_.tau_w + std::min(ode_state[State_::V_m], node.P_.V_peak) / node.P_.tau_w) - ode_state[State_::w] / node.P_.tau_w;
  f[State_::g_AMPA__X__AMPA] = (-ode_state[State_::g_AMPA__X__AMPA]) / node.P_.tau_decay_AMPA;
  f[State_::g_NMDA__X__NMDA] = ode_state[State_::g_NMDA__DOLLAR__X__NMDA] - ode_state[State_::g_NMDA__X__NMDA] / node.P_.tau_syn_rise_NMDA;
  f[State_::g_NMDA__DOLLAR__X__NMDA] = (-ode_state[State_::g_NMDA__DOLLAR__X__NMDA]) / node.P_.tau_syn_decay_NMDA;
  f[State_::g_GABA__X__GABA] = (-ode_state[State_::g_GABA__X__GABA]) / node.P_.tau_decay_GABA;
  return GSL_SUCCESS;
}

void adex_gamma_E_ml::update(nest::Time const & origin,const long from, const long to)
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function



  for ( long lag = from ; lag < to ; ++lag )
  {
    auto get_t = [origin, lag](){ return nest::Time( nest::Time::step( origin.get_steps() + lag + 1) ).get_ms(); };

    for (long i = 0; i < NUM_SPIKE_RECEPTORS; ++i)
    {
        get_spike_inputs_grid_sum_()[i] = get_spike_inputs_()[i].get_value(lag);
    }
    B_.I_stim_grid_sum_ = get_I_stim().get_value(lag);

    // NESTML generated code for the update block
  double g_AMPA__X__AMPA__tmp = V_.__P__g_AMPA__X__AMPA__g_AMPA__X__AMPA * S_.ode_state[State_::g_AMPA__X__AMPA];
  double g_NMDA__X__NMDA__tmp = V_.__P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA * S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA] + V_.__P__g_NMDA__X__NMDA__g_NMDA__X__NMDA * S_.ode_state[State_::g_NMDA__X__NMDA];
  double g_NMDA__DOLLAR__X__NMDA__tmp = V_.__P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA * S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA];
  double g_GABA__X__GABA__tmp = V_.__P__g_GABA__X__GABA__g_GABA__X__GABA * S_.ode_state[State_::g_GABA__X__GABA];
  double __t = 0;
  // numerical integration with adaptive step size control:
  // ------------------------------------------------------
  // gsl_odeiv_evolve_apply performs only a single numerical
  // integration step, starting from t and bounded by step;
  // the while-loop ensures integration over the whole simulation
  // step (0, step] if more than one integration step is needed due
  // to a small integration step size;
  // note that (t+IntegrationStep > step) leads to integration over
  // (t, step] and afterwards setting t to step, but it does not
  // enforce setting IntegrationStep to step-t; this is of advantage
  // for a consistent and efficient integration across subsequent
  // simulation intervals
  while ( __t < B_.__step )
  {

    const int status = gsl_odeiv_evolve_apply(B_.__e,
                                              B_.__c,
                                              B_.__s,
                                              &B_.__sys,              // system of ODE
                                              &__t,                   // from t
                                              B_.__step,              // to t <= step
                                              &B_.__integration_step, // integration step size
                                              S_.ode_state);          // neuronal state

    if ( status != GSL_SUCCESS )
    {
      throw nest::GSLSolverFailure( get_name(), status );
    }
  }
  /* replace analytically solvable variables with precisely integrated values  */
  S_.ode_state[State_::g_AMPA__X__AMPA] = g_AMPA__X__AMPA__tmp;
  S_.ode_state[State_::g_NMDA__X__NMDA] = g_NMDA__X__NMDA__tmp;
  S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA] = g_NMDA__DOLLAR__X__NMDA__tmp;
  S_.ode_state[State_::g_GABA__X__GABA] = g_GABA__X__GABA__tmp;
  S_.ode_state[State_::g_AMPA__X__AMPA] += ((0.001 * B_.spike_inputs_grid_sum_[AMPA - MIN_SPIKE_RECEPTOR])) / (1 / 1000.0);
  S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA] += ((0.001 * B_.spike_inputs_grid_sum_[NMDA - MIN_SPIKE_RECEPTOR])) * (V_.g_NMDA_const * (1 / P_.tau_syn_rise_NMDA - 1 / P_.tau_syn_decay_NMDA)) / (1 / 1000.0);
  S_.ode_state[State_::g_GABA__X__GABA] += ((0.001 * B_.spike_inputs_grid_sum_[GABA - MIN_SPIKE_RECEPTOR])) / (1 / 1000.0);
  if (S_.r > 0)
  {  
    S_.r -= 1;
    S_.ode_state[State_::V_m] = P_.V_reset;
  }
  else if (S_.ode_state[State_::V_m] >= P_.V_peak)
  {  
    S_.r = V_.RefractoryCounts;
    S_.ode_state[State_::V_m] = P_.V_reset;
    S_.ode_state[State_::w] += P_.b;
    set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
    nest::SpikeEvent se;
    nest::kernel().event_delivery_manager.send(*this, se, lag);
  }
    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void adex_gamma_E_ml::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


void adex_gamma_E_ml::handle(nest::SpikeEvent &e)
{
  assert(e.get_delay_steps() > 0);
  assert( e.get_rport() < B_.spike_inputs_.size() );

  double weight = e.get_weight();
  size_t nestml_buffer_idx = 0;
  if ( weight >= 0.0 )
  {
    nestml_buffer_idx = std::get<0>(rport_to_nestml_buffer_idx[e.get_rport()]);
  }
  else
  {
    nestml_buffer_idx = std::get<1>(rport_to_nestml_buffer_idx[e.get_rport()]);
    if ( nestml_buffer_idx == adex_gamma_E_ml::PORT_NOT_AVAILABLE )
    {
      nestml_buffer_idx = std::get<0>(rport_to_nestml_buffer_idx[e.get_rport()]);
    }
    weight = -weight;
  }
  B_.spike_inputs_[ nestml_buffer_idx - MIN_SPIKE_RECEPTOR ].add_value(
    e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
    weight * e.get_multiplicity() );
}

void adex_gamma_E_ml::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();     // we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_I_stim().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
}

