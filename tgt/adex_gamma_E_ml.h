
/**
 *  adex_gamma_E_ml.h
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
#ifndef ADEX_GAMMA_E_ML
#define ADEX_GAMMA_E_ML

#ifndef HAVE_LIBLTDL
#error "NEST was compiled without support for dynamic loading. Please install libltdl and recompile NEST."
#endif

// C++ includes:
#include <cmath>

#include "config.h"

#ifndef HAVE_GSL
#error "The GSL library is required for the Runge-Kutta solver."
#endif

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "dict_util.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{
namespace adex_gamma_E_ml_names
{
    const Name _V_m( "V_m" );
    const Name _w( "w" );
    const Name _r( "r" );
    const Name _g_AMPA__X__AMPA( "g_AMPA__X__AMPA" );
    const Name _g_NMDA__X__NMDA( "g_NMDA__X__NMDA" );
    const Name _g_NMDA__DOLLAR__X__NMDA( "g_NMDA__DOLLAR__X__NMDA" );
    const Name _g_GABA__X__GABA( "g_GABA__X__GABA" );
    const Name _I_AMPA( "I_AMPA" );
    const Name _I_GABA( "I_GABA" );
    const Name _I_NMDA( "I_NMDA" );
    const Name _C_m( "C_m" );
    const Name _t_ref( "t_ref" );
    const Name _V_reset( "V_reset" );
    const Name _g_L( "g_L" );
    const Name _E_L( "E_L" );
    const Name _a( "a" );
    const Name _b( "b" );
    const Name _Delta_T( "Delta_T" );
    const Name _tau_w( "tau_w" );
    const Name _V_th( "V_th" );
    const Name _V_peak( "V_peak" );
    const Name _E_e( "E_e" );
    const Name _E_i( "E_i" );
    const Name _tau_decay_AMPA( "tau_decay_AMPA" );
    const Name _tau_decay_GABA( "tau_decay_GABA" );
    const Name _Mg( "Mg" );
    const Name _tau_syn_rise_NMDA( "tau_syn_rise_NMDA" );
    const Name _tau_syn_decay_NMDA( "tau_syn_decay_NMDA" );
    const Name _I_e( "I_e" );
}
}



/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
**/
extern "C" inline int adex_gamma_E_ml_dynamics( double, const double y[], double f[], void* pnode );


#include "nest_time.h"
  typedef nest::port nest_port_t;
  typedef nest::rport nest_rport_t;

/* BeginDocumentation
  Name: adex_gamma_E_ml

  Description:

    

  Parameters:
  The following parameters can be set in the status dictionary.
C_m [pF]  membrane parameters
 Membrane Capacitance
t_ref [ms]  Refractory period
V_reset [mV]  Reset Potential
g_L [nS]  Leak Conductance
E_L [mV]  Leak reversal Potential (aka resting potential)
a [nS]  spike adaptation parameters
 Subthreshold adaptation
b [pA]  Spike-triggered adaptation
Delta_T [mV]  Slope factor
tau_w [ms]  Adaptation time constant
V_th [mV]  Threshold Potential
V_peak [mV]  Spike detection threshold
tau_decay_AMPA [ms]  Time decay constants
tau_syn_rise_NMDA [ms]  Synaptic time constant excitatory synapse
I_e [pA]  constant external input current


  Dynamic state variables:
V_m [mV]  Membrane potential
w [pA]  Spike-adaptation current
r [integer]  Counts number of tick during the refractory period


  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class adex_gamma_E_ml : public nest::ArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  adex_gamma_E_ml();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook() (or calibrate() in NEST 3.3 and older).
  **/
  adex_gamma_E_ml(const adex_gamma_E_ml &);

  /**
   * Destructor.
  **/
  ~adex_gamma_E_ml() override;

  // -------------------------------------------------------------------------
  //   Import sets of overloaded virtual functions.
  //   See: Technical Issues / Virtual Functions: Overriding, Overloading,
  //        and Hiding
  // -------------------------------------------------------------------------

  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
   * Used to validate that we can send nest::SpikeEvent to desired target:port.
  **/
  nest_port_t send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool) override;


  // -------------------------------------------------------------------------
  //   Functions handling incoming events.
  //   We tell nest that we can handle incoming events of various types by
  //   defining handle() for the given event.
  // -------------------------------------------------------------------------


  void handle(nest::SpikeEvent &) override;        //! accept spikes
  void handle(nest::CurrentEvent &) override;      //! accept input current

  void handle(nest::DataLoggingRequest &) override;//! allow recording with multimeter
  nest_port_t handles_test_event(nest::SpikeEvent&, nest_port_t) override;
  nest_port_t handles_test_event(nest::CurrentEvent&, nest_port_t) override;
  nest_port_t handles_test_event(nest::DataLoggingRequest&, nest_port_t) override;

  // -------------------------------------------------------------------------
  //   Functions for getting/setting parameters and state values.
  // -------------------------------------------------------------------------

  void get_status(DictionaryDatum &) const override;
  void set_status(const DictionaryDatum &) override;


  // -------------------------------------------------------------------------
  //   Getters/setters for state block
  // -------------------------------------------------------------------------

  inline double get_V_m() const
  {
    return S_.ode_state[State_::V_m];
  }

  inline void set_V_m(const double __v)
  {
    S_.ode_state[State_::V_m] = __v;
  }

  inline double get_w() const
  {
    return S_.ode_state[State_::w];
  }

  inline void set_w(const double __v)
  {
    S_.ode_state[State_::w] = __v;
  }

  inline long get_r() const
  {
    return S_.r;
  }

  inline void set_r(const long __v)
  {
    S_.r = __v;
  }

  inline double get_g_AMPA__X__AMPA() const
  {
    return S_.ode_state[State_::g_AMPA__X__AMPA];
  }

  inline void set_g_AMPA__X__AMPA(const double __v)
  {
    S_.ode_state[State_::g_AMPA__X__AMPA] = __v;
  }

  inline double get_g_NMDA__X__NMDA() const
  {
    return S_.ode_state[State_::g_NMDA__X__NMDA];
  }

  inline void set_g_NMDA__X__NMDA(const double __v)
  {
    S_.ode_state[State_::g_NMDA__X__NMDA] = __v;
  }

  inline double get_g_NMDA__DOLLAR__X__NMDA() const
  {
    return S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA];
  }

  inline void set_g_NMDA__DOLLAR__X__NMDA(const double __v)
  {
    S_.ode_state[State_::g_NMDA__DOLLAR__X__NMDA] = __v;
  }

  inline double get_g_GABA__X__GABA() const
  {
    return S_.ode_state[State_::g_GABA__X__GABA];
  }

  inline void set_g_GABA__X__GABA(const double __v)
  {
    S_.ode_state[State_::g_GABA__X__GABA] = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------

  inline double get_C_m() const
  {
    return P_.C_m;
  }

  inline void set_C_m(const double __v)
  {
    P_.C_m = __v;
  }

  inline double get_t_ref() const
  {
    return P_.t_ref;
  }

  inline void set_t_ref(const double __v)
  {
    P_.t_ref = __v;
  }

  inline double get_V_reset() const
  {
    return P_.V_reset;
  }

  inline void set_V_reset(const double __v)
  {
    P_.V_reset = __v;
  }

  inline double get_g_L() const
  {
    return P_.g_L;
  }

  inline void set_g_L(const double __v)
  {
    P_.g_L = __v;
  }

  inline double get_E_L() const
  {
    return P_.E_L;
  }

  inline void set_E_L(const double __v)
  {
    P_.E_L = __v;
  }

  inline double get_a() const
  {
    return P_.a;
  }

  inline void set_a(const double __v)
  {
    P_.a = __v;
  }

  inline double get_b() const
  {
    return P_.b;
  }

  inline void set_b(const double __v)
  {
    P_.b = __v;
  }

  inline double get_Delta_T() const
  {
    return P_.Delta_T;
  }

  inline void set_Delta_T(const double __v)
  {
    P_.Delta_T = __v;
  }

  inline double get_tau_w() const
  {
    return P_.tau_w;
  }

  inline void set_tau_w(const double __v)
  {
    P_.tau_w = __v;
  }

  inline double get_V_th() const
  {
    return P_.V_th;
  }

  inline void set_V_th(const double __v)
  {
    P_.V_th = __v;
  }

  inline double get_V_peak() const
  {
    return P_.V_peak;
  }

  inline void set_V_peak(const double __v)
  {
    P_.V_peak = __v;
  }

  inline double get_E_e() const
  {
    return P_.E_e;
  }

  inline void set_E_e(const double __v)
  {
    P_.E_e = __v;
  }

  inline double get_E_i() const
  {
    return P_.E_i;
  }

  inline void set_E_i(const double __v)
  {
    P_.E_i = __v;
  }

  inline double get_tau_decay_AMPA() const
  {
    return P_.tau_decay_AMPA;
  }

  inline void set_tau_decay_AMPA(const double __v)
  {
    P_.tau_decay_AMPA = __v;
  }

  inline double get_tau_decay_GABA() const
  {
    return P_.tau_decay_GABA;
  }

  inline void set_tau_decay_GABA(const double __v)
  {
    P_.tau_decay_GABA = __v;
  }

  inline double get_Mg() const
  {
    return P_.Mg;
  }

  inline void set_Mg(const double __v)
  {
    P_.Mg = __v;
  }

  inline double get_tau_syn_rise_NMDA() const
  {
    return P_.tau_syn_rise_NMDA;
  }

  inline void set_tau_syn_rise_NMDA(const double __v)
  {
    P_.tau_syn_rise_NMDA = __v;
  }

  inline double get_tau_syn_decay_NMDA() const
  {
    return P_.tau_syn_decay_NMDA;
  }

  inline void set_tau_syn_decay_NMDA(const double __v)
  {
    P_.tau_syn_decay_NMDA = __v;
  }

  inline double get_I_e() const
  {
    return P_.I_e;
  }

  inline void set_I_e(const double __v)
  {
    P_.I_e = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for internals
  // -------------------------------------------------------------------------

  inline double get_t_peak_NMDA() const
  {
    return V_.t_peak_NMDA;
  }

  inline void set_t_peak_NMDA(const double __v)
  {
    V_.t_peak_NMDA = __v;
  }
  inline double get___h() const
  {
    return V_.__h;
  }

  inline void set___h(const double __v)
  {
    V_.__h = __v;
  }
  inline double get_g_NMDA_const() const
  {
    return V_.g_NMDA_const;
  }

  inline void set_g_NMDA_const(const double __v)
  {
    V_.g_NMDA_const = __v;
  }
  inline long get_RefractoryCounts() const
  {
    return V_.RefractoryCounts;
  }

  inline void set_RefractoryCounts(const long __v)
  {
    V_.RefractoryCounts = __v;
  }
  inline double get___P__g_AMPA__X__AMPA__g_AMPA__X__AMPA() const
  {
    return V_.__P__g_AMPA__X__AMPA__g_AMPA__X__AMPA;
  }

  inline void set___P__g_AMPA__X__AMPA__g_AMPA__X__AMPA(const double __v)
  {
    V_.__P__g_AMPA__X__AMPA__g_AMPA__X__AMPA = __v;
  }
  inline double get___P__g_NMDA__X__NMDA__g_NMDA__X__NMDA() const
  {
    return V_.__P__g_NMDA__X__NMDA__g_NMDA__X__NMDA;
  }

  inline void set___P__g_NMDA__X__NMDA__g_NMDA__X__NMDA(const double __v)
  {
    V_.__P__g_NMDA__X__NMDA__g_NMDA__X__NMDA = __v;
  }
  inline double get___P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA() const
  {
    return V_.__P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA;
  }

  inline void set___P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA(const double __v)
  {
    V_.__P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA = __v;
  }
  inline double get___P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA() const
  {
    return V_.__P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA;
  }

  inline void set___P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA(const double __v)
  {
    V_.__P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA = __v;
  }
  inline double get___P__g_GABA__X__GABA__g_GABA__X__GABA() const
  {
    return V_.__P__g_GABA__X__GABA__g_GABA__X__GABA;
  }

  inline void set___P__g_GABA__X__GABA__g_GABA__X__GABA(const double __v)
  {
    V_.__P__g_GABA__X__GABA__g_GABA__X__GABA = __v;
  }


  // -------------------------------------------------------------------------
  //   Initialization functions
  // -------------------------------------------------------------------------
  void calibrate_time( const nest::TimeConverter& tc ) override;

protected:

private:
  void recompute_internal_variables(bool exclude_timestep=false);

private:
/**
   * Synapse types to connect to
   * @note Excluded lower and upper bounds are defined as MIN_, MAX_.
   *       Excluding port 0 avoids accidental connections.
  **/
  static const nest_port_t MIN_SPIKE_RECEPTOR = 1;
  static const nest_port_t PORT_NOT_AVAILABLE = -1;

  enum SynapseTypes
  {
    AMPA = 1,
    NMDA = 2,
    GABA = 3,
    MAX_SPIKE_RECEPTOR = 4
  };

  static const size_t NUM_SPIKE_RECEPTORS = MAX_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR;

static std::vector< std::tuple< int, int > > rport_to_nestml_buffer_idx;

  /**
   * Reset state of neuron.
  **/

  void init_state_internal_();

  /**
   * Reset internal buffers of neuron.
  **/
  void init_buffers_() override;

  /**
   * Initialize auxiliary quantities, leave parameters and state untouched.
  **/
  void pre_run_hook() override;

  /**
   * Take neuron through given time interval
  **/
  void update(nest::Time const &, const long, const long) override;

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap<adex_gamma_E_ml>;
  friend class nest::UniversalDataLogger<adex_gamma_E_ml>;

  /**
   * Free parameters of the neuron.
   *

 parameters done
   *
   * These are the parameters that can be set by the user through @c `node.set()`.
   * They are initialized from the model prototype when the node is created.
   * Parameters do not change during calls to @c update() and are not reset by
   * @c ResetNetwork.
   *
   * @note Parameters_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If Parameters_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If Parameters_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct Parameters_
  {    
    //!  membrane parameters
    //!  Membrane Capacitance
    double C_m;
    //!  Refractory period
    double t_ref;
    //!  Reset Potential
    double V_reset;
    //!  Leak Conductance
    double g_L;
    //!  Leak reversal Potential (aka resting potential)
    double E_L;
    //!  spike adaptation parameters
    //!  Subthreshold adaptation
    double a;
    //!  Spike-triggered adaptation
    double b;
    //!  Slope factor
    double Delta_T;
    //!  Adaptation time constant
    double tau_w;
    //!  Threshold Potential
    double V_th;
    //!  Spike detection threshold
    double V_peak;
    double E_e;
    double E_i;
    //!  Time decay constants
    double tau_decay_AMPA;
    double tau_decay_GABA;
    double Mg;
    //!  Synaptic time constant excitatory synapse
    double tau_syn_rise_NMDA;
    double tau_syn_decay_NMDA;
    //!  constant external input current
    double I_e;

    double __gsl_error_tol;

    /**
     * Initialize parameters to their default values.
    **/
    Parameters_();
  };

  /**
   * Dynamic state of the neuron.
   *
   *
   *
   * These are the state variables that are advanced in time by calls to
   * @c update(). In many models, some or all of them can be set by the user
   * through @c `node.set()`. The state variables are initialized from the model
   * prototype when the node is created. State variables are reset by @c ResetNetwork.
   *
   * @note State_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If State_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If State_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct State_
  {

    // non-ODE state variables
//!  Counts number of tick during the refractory period
long r;
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_m,
      w,
      g_AMPA__X__AMPA,
      g_NMDA__X__NMDA,
      g_NMDA__DOLLAR__X__NMDA,
      g_GABA__X__GABA,
      // moved state variables from synapse (numeric)
      // moved state variables from synapse (analytic)
      // final entry to easily get the vector size
      STATE_VEC_SIZE
    };

    //! state vector, must be C-array for GSL solver
    double ode_state[STATE_VEC_SIZE];

    State_();
  };

  struct DelayedVariables_
  {
  };

  /**
   * Internal variables of the neuron.
   *
   *
   *
   * These variables must be initialized by @c pre_run_hook (or calibrate in NEST 3.3 and older), which is called before
   * the first call to @c update() upon each call to @c Simulate.
   * @node Variables_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c pre_run_hook() (or calibrate() in NEST 3.3 and older). If Variables_ has members that
   *       cannot destroy themselves, Variables_ will need a destructor.
  **/
  struct Variables_
  {
    double t_peak_NMDA;
    double __h;
    double g_NMDA_const;
    long RefractoryCounts;
    double __P__g_AMPA__X__AMPA__g_AMPA__X__AMPA;
    double __P__g_NMDA__X__NMDA__g_NMDA__X__NMDA;
    double __P__g_NMDA__X__NMDA__g_NMDA__DOLLAR__X__NMDA;
    double __P__g_NMDA__DOLLAR__X__NMDA__g_NMDA__DOLLAR__X__NMDA;
    double __P__g_GABA__X__GABA__g_GABA__X__GABA;
  };

  /**
   * Buffers of the neuron.
   * Usually buffers for incoming spikes and data logged for analog recorders.
   * Buffers must be initialized by @c init_buffers_(), which is called before
   * @c pre_run_hook() (or calibrate() in NEST 3.3 and older) on the first call to @c Simulate after the start of NEST,
   * ResetKernel or ResetNetwork.
   * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
   *       cannot destroy themselves, Buffers_ will need a destructor.
  **/
  struct Buffers_
  {
    Buffers_(adex_gamma_E_ml &);
    Buffers_(const Buffers_ &, adex_gamma_E_ml &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<adex_gamma_E_ml> logger_;

    // -----------------------------------------------------------------------
    //   Buffers and sums of incoming spikes/currents per timestep
    // -----------------------------------------------------------------------
    // Buffer containing the incoming spikes
    

inline std::vector< nest::RingBuffer >& get_spike_inputs_()
{
    return spike_inputs_;
}
std::vector< nest::RingBuffer > spike_inputs_;

    // Buffer containing the sum of all the incoming spikes
    

inline std::vector< double >& get_spike_inputs_grid_sum_()
{
    return spike_inputs_grid_sum_;
}
std::vector< double > spike_inputs_grid_sum_;

nest::RingBuffer
 I_stim;   //!< Buffer for input (type: pA)    
    inline nest::RingBuffer& get_I_stim() {
        return I_stim;
    }

double I_stim_grid_sum_;

    // -----------------------------------------------------------------------
    //   GSL ODE solver data structures
    // -----------------------------------------------------------------------

    gsl_odeiv_step* __s;    //!< stepping function
    gsl_odeiv_control* __c; //!< adaptive stepsize control function
    gsl_odeiv_evolve* __e;  //!< evolution function
    gsl_odeiv_system __sys; //!< struct describing system

    // __integration_step should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double __step;             //!< step size in ms
    double __integration_step; //!< current integration time step, updated by GSL
  };

  // -------------------------------------------------------------------------
  //   Getters/setters for inline expressions
  // -------------------------------------------------------------------------
  inline double get_V_bounded() const
  {
    return std::min(S_.ode_state[State_::V_m], P_.V_peak);
  }

  inline double get_B() const
  {
    return 1 / (1 + std::exp((-0.062) * (std::min(S_.ode_state[State_::V_m], P_.V_peak)) * 1 / 1.0) * (P_.Mg / 3.57));
  }

  inline double get_I_AMPA() const
  {
    return S_.ode_state[State_::g_AMPA__X__AMPA] * 1.0 * ((std::min(S_.ode_state[State_::V_m], P_.V_peak)) - P_.E_e);
  }

  inline double get_I_GABA() const
  {
    return S_.ode_state[State_::g_GABA__X__GABA] * 1.0 * ((std::min(S_.ode_state[State_::V_m], P_.V_peak)) - P_.E_i);
  }

  inline double get_I_NMDA() const
  {
    return S_.ode_state[State_::g_NMDA__X__NMDA] * 1.0 * ((std::min(S_.ode_state[State_::V_m], P_.V_peak)) - P_.E_e) * (1 / (1 + std::exp((-0.062) * (std::min(S_.ode_state[State_::V_m], P_.V_peak)) * 1 / 1.0) * (P_.Mg / 3.57)));
  }

  inline double get_I_syn() const
  {
    return (S_.ode_state[State_::g_AMPA__X__AMPA] * 1.0 * ((std::min(S_.ode_state[State_::V_m], P_.V_peak)) - P_.E_e)) + (S_.ode_state[State_::g_NMDA__X__NMDA] * 1.0 * ((std::min(S_.ode_state[State_::V_m], P_.V_peak)) - P_.E_e) * (1 / (1 + std::exp((-0.062) * (std::min(S_.ode_state[State_::V_m], P_.V_peak)) * 1 / 1.0) * (P_.Mg / 3.57)))) + (S_.ode_state[State_::g_GABA__X__GABA] * 1.0 * ((std::min(S_.ode_state[State_::V_m], P_.V_peak)) - P_.E_i));
  }



  // -------------------------------------------------------------------------
  //   Getters/setters for input buffers
  // -------------------------------------------------------------------------

  // Buffer containing the incoming spikes
  

inline std::vector< nest::RingBuffer >& get_spike_inputs_()
{
    return B_.get_spike_inputs_();
}

  

inline std::vector< double >& get_spike_inputs_grid_sum_()
{
    return B_.get_spike_inputs_grid_sum_();
}
  
inline nest::RingBuffer& get_I_stim() {
    return B_.get_I_stim();
}

  // -------------------------------------------------------------------------
  //   Member variables of neuron model.
  //   Each model neuron should have precisely the following four data members,
  //   which are one instance each of the parameters, state, buffers and variables
  //   structures. Experience indicates that the state and variables member should
  //   be next to each other to achieve good efficiency (caching).
  //   Note: Devices require one additional data member, an instance of the
  //   ``Device`` child class they belong to.
  // -------------------------------------------------------------------------


  Parameters_       P_;        //!< Free parameters.
  State_            S_;        //!< Dynamic state.
  DelayedVariables_ DV_;       //!< Delayed state variables.
  Variables_        V_;        //!< Internal Variables
  Buffers_          B_;        //!< Buffers.

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap<adex_gamma_E_ml> recordablesMap_;
  friend int adex_gamma_E_ml_dynamics( double, const double y[], double f[], void* pnode );


}; /* neuron adex_gamma_E_ml */

inline nest_port_t adex_gamma_E_ml::send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest_port_t adex_gamma_E_ml::handles_test_event(nest::SpikeEvent&, nest_port_t receptor_type)
{
    assert( B_.spike_inputs_.size() == NUM_SPIKE_RECEPTORS );
    if ( receptor_type < MIN_SPIKE_RECEPTOR or receptor_type >= MAX_SPIKE_RECEPTOR )
    {
      throw nest::UnknownReceptorType( receptor_type, get_name() );
    }
    return receptor_type - MIN_SPIKE_RECEPTOR;
}

inline nest_port_t adex_gamma_E_ml::handles_test_event(nest::CurrentEvent&, nest_port_t receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }
  return 0;
}

inline nest_port_t adex_gamma_E_ml::handles_test_event(nest::DataLoggingRequest& dlr, nest_port_t receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline void adex_gamma_E_ml::get_status(DictionaryDatum &__d) const
{
  // parameters
  def<double>(__d, nest::adex_gamma_E_ml_names::_C_m, get_C_m());
  def<double>(__d, nest::adex_gamma_E_ml_names::_t_ref, get_t_ref());
  def<double>(__d, nest::adex_gamma_E_ml_names::_V_reset, get_V_reset());
  def<double>(__d, nest::adex_gamma_E_ml_names::_g_L, get_g_L());
  def<double>(__d, nest::adex_gamma_E_ml_names::_E_L, get_E_L());
  def<double>(__d, nest::adex_gamma_E_ml_names::_a, get_a());
  def<double>(__d, nest::adex_gamma_E_ml_names::_b, get_b());
  def<double>(__d, nest::adex_gamma_E_ml_names::_Delta_T, get_Delta_T());
  def<double>(__d, nest::adex_gamma_E_ml_names::_tau_w, get_tau_w());
  def<double>(__d, nest::adex_gamma_E_ml_names::_V_th, get_V_th());
  def<double>(__d, nest::adex_gamma_E_ml_names::_V_peak, get_V_peak());
  def<double>(__d, nest::adex_gamma_E_ml_names::_E_e, get_E_e());
  def<double>(__d, nest::adex_gamma_E_ml_names::_E_i, get_E_i());
  def<double>(__d, nest::adex_gamma_E_ml_names::_tau_decay_AMPA, get_tau_decay_AMPA());
  def<double>(__d, nest::adex_gamma_E_ml_names::_tau_decay_GABA, get_tau_decay_GABA());
  def<double>(__d, nest::adex_gamma_E_ml_names::_Mg, get_Mg());
  def<double>(__d, nest::adex_gamma_E_ml_names::_tau_syn_rise_NMDA, get_tau_syn_rise_NMDA());
  def<double>(__d, nest::adex_gamma_E_ml_names::_tau_syn_decay_NMDA, get_tau_syn_decay_NMDA());
  def<double>(__d, nest::adex_gamma_E_ml_names::_I_e, get_I_e());

  // initial values for state variables in ODE or kernel
  def<double>(__d, nest::adex_gamma_E_ml_names::_V_m, get_V_m());
  def<double>(__d, nest::adex_gamma_E_ml_names::_w, get_w());
  def<long>(__d, nest::adex_gamma_E_ml_names::_r, get_r());
  def<double>(__d, nest::adex_gamma_E_ml_names::_g_AMPA__X__AMPA, get_g_AMPA__X__AMPA());
  def<double>(__d, nest::adex_gamma_E_ml_names::_g_NMDA__X__NMDA, get_g_NMDA__X__NMDA());
  def<double>(__d, nest::adex_gamma_E_ml_names::_g_NMDA__DOLLAR__X__NMDA, get_g_NMDA__DOLLAR__X__NMDA());
  def<double>(__d, nest::adex_gamma_E_ml_names::_g_GABA__X__GABA, get_g_GABA__X__GABA());

  ArchivingNode::get_status( __d );
  DictionaryDatum __receptor_type = new Dictionary();
    ( *__receptor_type )[ "AMPA" ] = 1;
    ( *__receptor_type )[ "NMDA" ] = 2;
    ( *__receptor_type )[ "GABA" ] = 3;
    ( *__d )[ "receptor_types" ] = __receptor_type;

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

inline void adex_gamma_E_ml::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_C_m = get_C_m();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_C_m, tmp_C_m, this);
  double tmp_t_ref = get_t_ref();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_t_ref, tmp_t_ref, this);
  double tmp_V_reset = get_V_reset();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_V_reset, tmp_V_reset, this);
  double tmp_g_L = get_g_L();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_g_L, tmp_g_L, this);
  double tmp_E_L = get_E_L();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_E_L, tmp_E_L, this);
  double tmp_a = get_a();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_a, tmp_a, this);
  double tmp_b = get_b();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_b, tmp_b, this);
  double tmp_Delta_T = get_Delta_T();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_Delta_T, tmp_Delta_T, this);
  double tmp_tau_w = get_tau_w();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_tau_w, tmp_tau_w, this);
  double tmp_V_th = get_V_th();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_V_th, tmp_V_th, this);
  double tmp_V_peak = get_V_peak();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_V_peak, tmp_V_peak, this);
  double tmp_E_e = get_E_e();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_E_e, tmp_E_e, this);
  double tmp_E_i = get_E_i();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_E_i, tmp_E_i, this);
  double tmp_tau_decay_AMPA = get_tau_decay_AMPA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_tau_decay_AMPA, tmp_tau_decay_AMPA, this);
  double tmp_tau_decay_GABA = get_tau_decay_GABA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_tau_decay_GABA, tmp_tau_decay_GABA, this);
  double tmp_Mg = get_Mg();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_Mg, tmp_Mg, this);
  double tmp_tau_syn_rise_NMDA = get_tau_syn_rise_NMDA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_tau_syn_rise_NMDA, tmp_tau_syn_rise_NMDA, this);
  double tmp_tau_syn_decay_NMDA = get_tau_syn_decay_NMDA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_tau_syn_decay_NMDA, tmp_tau_syn_decay_NMDA, this);
  double tmp_I_e = get_I_e();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_I_e, tmp_I_e, this);

  // initial values for state variables in ODE or kernel
  double tmp_V_m = get_V_m();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_V_m, tmp_V_m, this);
  double tmp_w = get_w();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_w, tmp_w, this);
  long tmp_r = get_r();
  nest::updateValueParam<long>(__d, nest::adex_gamma_E_ml_names::_r, tmp_r, this);
  double tmp_g_AMPA__X__AMPA = get_g_AMPA__X__AMPA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_g_AMPA__X__AMPA, tmp_g_AMPA__X__AMPA, this);
  double tmp_g_NMDA__X__NMDA = get_g_NMDA__X__NMDA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_g_NMDA__X__NMDA, tmp_g_NMDA__X__NMDA, this);
  double tmp_g_NMDA__DOLLAR__X__NMDA = get_g_NMDA__DOLLAR__X__NMDA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_g_NMDA__DOLLAR__X__NMDA, tmp_g_NMDA__DOLLAR__X__NMDA, this);
  double tmp_g_GABA__X__GABA = get_g_GABA__X__GABA();
  nest::updateValueParam<double>(__d, nest::adex_gamma_E_ml_names::_g_GABA__X__GABA, tmp_g_GABA__X__GABA, this);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_C_m(tmp_C_m);
  set_t_ref(tmp_t_ref);
  set_V_reset(tmp_V_reset);
  set_g_L(tmp_g_L);
  set_E_L(tmp_E_L);
  set_a(tmp_a);
  set_b(tmp_b);
  set_Delta_T(tmp_Delta_T);
  set_tau_w(tmp_tau_w);
  set_V_th(tmp_V_th);
  set_V_peak(tmp_V_peak);
  set_E_e(tmp_E_e);
  set_E_i(tmp_E_i);
  set_tau_decay_AMPA(tmp_tau_decay_AMPA);
  set_tau_decay_GABA(tmp_tau_decay_GABA);
  set_Mg(tmp_Mg);
  set_tau_syn_rise_NMDA(tmp_tau_syn_rise_NMDA);
  set_tau_syn_decay_NMDA(tmp_tau_syn_decay_NMDA);
  set_I_e(tmp_I_e);
  set_V_m(tmp_V_m);
  set_w(tmp_w);
  set_r(tmp_r);
  set_g_AMPA__X__AMPA(tmp_g_AMPA__X__AMPA);
  set_g_NMDA__X__NMDA(tmp_g_NMDA__X__NMDA);
  set_g_NMDA__DOLLAR__X__NMDA(tmp_g_NMDA__DOLLAR__X__NMDA);
  set_g_GABA__X__GABA(tmp_g_GABA__X__GABA);




  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }

  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
};



#endif /* #ifndef ADEX_GAMMA_E_ML */
