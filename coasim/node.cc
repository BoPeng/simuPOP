/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "node.hh"

#ifndef CORE__COMPILE_OPTIONS_HH
# include "compile_options.hh"
#endif

#ifndef CORE__DIST_FUNCTIONS_HH_INCLUDED
# include "dist_funcs.hh"
#endif
#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif

#ifndef SSTREAM_INCLUDED
# include <sstream>
# define SSTREAM_INCLUDED
#endif
#ifndef STRING_INCLUDED
# include <string>
# define STRING_INCLUDED
#endif

#ifndef CASSERT_INCLUDED
# include <cassert>
# define CASSERT_INCLUDED
#endif
#if EXPENSIVE_ASSERTS
# ifndef FSTREAM_INCLUDED
#  include <fstream>
#  define FSTREAM_INCLUDED
# endif
#endif

#ifndef CSTDLIB_INCLUDED
# include <cstdlib>
# define CSTDLIB_INCLUDED
#endif

namespace core
{

  void Node::initialize_marker(unsigned int idx, const Marker &m)
  {
    if (m_states.size() <= idx)
      throw std::out_of_range("marker index out of range");
    m_states[idx] = m.default_value();
  }

  double LeafNode::surface_at_point(double point) const
    throw(std::out_of_range)
  {
    if (point < 0 or 1.0 <= point)
      throw std::out_of_range("Point out of range [0,1).");
    return 0.0;
  }

  void LeafNode::print_tree_at_point(std::ostream &os, double point,
    double edge_length,
    bool print_edge) const
    throw(std::out_of_range)
  {
    if (point < 0 or 1.0 <= point)
      throw std::out_of_range("Point out of range [0,1).");
    os << '\'' << m_id << '\'';
    if (print_edge) os << " : " << edge_length;
  }

  void LeafNode::mutate_marker(unsigned int idx, Mutator &m)
  {
    // no mutations out of leaf
  }

  double CoalescentNode::surface_at_point(double point) const
    throw(std::out_of_range)
  {
    // NB! don't check if this node contains it -- it could be
    // retired, if that is the case it's children will contain it.
    // if this node does not contain it, neither of it's children
    // will, so it still works out (althoug we call recursively at
    // little more than strictly necessary)
    double surface = 0.0;
    if (m_left->intervals().contains_point(point))
    {
      surface += m_left->surface_at_point(point);
      surface += time() - m_left->time();
    }
    if (m_right->intervals().contains_point(point))
    {
      surface += m_right->surface_at_point(point);
      surface += time() - m_right->time();
    }
    return surface;
  }

  void CoalescentNode::print_tree_at_point(std::ostream &os, double point,
    double edge_length,
    bool print_edge) const
    throw(std::out_of_range)
  {
    double left_dist  = time() - m_left->time();
    double right_dist = time() - m_right->time();

    if (m_left->intervals().contains_point(point)
      and m_right->intervals().contains_point(point))
    {
      os << '(';
      m_left->print_tree_at_point(os, point, left_dist, true);
      os << ',';
      m_right->print_tree_at_point(os, point, right_dist, true);
      os << ')';
      if (print_edge) os << " : " << edge_length;

    }
    else
    {
      if (m_left->intervals().contains_point(point))
        m_left->print_tree_at_point(os, point,
          edge_length+left_dist,
          print_edge);
      if (m_right->intervals().contains_point(point))
        m_right->print_tree_at_point(os, point,
          edge_length+right_dist,
          print_edge);
    }
  }

  void CoalescentNode::mutate_marker(unsigned int idx, Mutator &m)
  {
    if (! (idx < no_states()) )
      throw std::out_of_range("marker index out of range");

    double point = m_conf.position(idx);

    if (m_left->intervals().contains_point(point))
    {
      set_state(m_left, idx, m.mutate(*this,*m_left,state(idx)));
      m_left->mutate_marker(idx,m);
    }

    if (m_right->intervals().contains_point(point))
    {
      set_state(m_right, idx, m.mutate(*this,*m_right,state(idx)));
      m_right->mutate_marker(idx,m);
    }
  }

  double RecombinationNode::surface_at_point(double point) const
    throw(std::out_of_range)
  {
    // no need to check here, it is the parents responsibility to
    //check that
    //if (! intervals().contains_point(point)) return 0.0;

    double surface = 0.0;
    if (m_child->intervals().contains_point(point))
    {
      surface += m_child->surface_at_point(point);
      surface += time() - m_child->time();
    }
    return surface;
  }

  void RecombinationNode::print_tree_at_point(std::ostream &os, double point,
    double edge_length,
    bool print_edge) const
    throw(std::out_of_range)
  {
    double d = time() - m_child->time();
    m_child->print_tree_at_point(os, point, edge_length+d, print_edge);
  }

  void RecombinationNode::mutate_marker(unsigned int idx, Mutator &m)
  {
    if (! (idx < no_states()) )
      throw std::out_of_range("marker index out of range");

    set_state(m_child, idx, m.mutate(*this,*m_child,state(idx)));
    m_child->mutate_marker(idx,m);
  }

  double GeneConversionNode::surface_at_point(double point) const
    throw(std::out_of_range)
  {
    // no need to check here, it is the parents responsibility to
    //check that
    //if (! intervals().contains_point(point)) return 0.0;

    double surface = 0.0;
    if (m_child->intervals().contains_point(point))
    {
      surface += m_child->surface_at_point(point);
      surface += time() - m_child->time();
    }
    return surface;
  }

  void GeneConversionNode::print_tree_at_point(std::ostream &os,
    double point,
    double edge_length,
    bool print_edge) const
    throw(std::out_of_range)
  {
    double d = time() - m_child->time();
    m_child->print_tree_at_point(os, point, edge_length+d, print_edge);
  }

  void GeneConversionNode::mutate_marker(unsigned int idx, Mutator &m)
  {
    if (! (idx < no_states()) )
      throw std::out_of_range("marker index out of range");

    set_state(m_child, idx, m.mutate(*this,*m_child,state(idx)));
    m_child->mutate_marker(idx,m);
  }

  ARG::~ARG()
  {
    std::vector<Node*>::iterator itr;
    for (itr = m_leaf_pool.begin(); itr != m_leaf_pool.end(); ++itr)
      delete *itr;
    for (itr = m_node_pool.begin(); itr != m_node_pool.end(); ++itr)
      delete *itr;
  }

  LeafNode *ARG::leaf() throw()
  {
    LeafNode *n = new LeafNode(m_conf, m_no_leaves);
    n->m_intervals.add(0.0,1.0,1);
    m_leaf_pool.push_back(n);
    ++m_no_leaves;

    return n;
  }

  static inline bool contains_marker(const Configuration &conf,
    const Interval &i)
  {
    for (int m = 0; m < conf.no_markers(); ++m)
      if (i.contains_point(conf.position(m))) return true;
    return false;
  }

  static Intervals filter_contains_marker(const Intervals &intervals,
    const Configuration &conf)
  {
    Intervals result;
    for (int i = 0; i < intervals.size(); ++i)
    {
      const Interval &it = intervals.interval(i);
      if (contains_marker(conf,it)) result.add(it);
    }
    return result;
  }

  CoalescentNode *ARG::coalescence(double time, Node *left, Node *right)
    throw(null_child)
  {
    if (left == 0 or right == 0) throw null_child();

    // sort in retired and non-retired intervals
    Intervals retired;
    Intervals non_retired;
    Intervals merged = left->intervals() | right->intervals();
    if (!m_keep_empty) merged = filter_contains_marker(merged, m_conf);

#if 0
    std::cout << "coalescence -- left: " << left->intervals() << std::endl;
    std::cout << "coalescence -- right: " << right->intervals() << std::endl;
    std::cout << "coalescence -- merged: " << merged << std::endl;
#endif

    for (int i = 0; i < merged.size(); ++i)
    {
      if (merged.interval(i).leaf_contacts() < m_no_leaves)
        non_retired.add(merged.interval(i));
      else if (merged.interval(i).leaf_contacts() == m_no_leaves)
        retired.add(merged.interval(i));
      else
        assert(false);
    }

    CoalescentNode *n = new CoalescentNode(m_conf,time,left,right,
      non_retired, retired);

    m_node_pool.push_back(n);

    for (int i = 0; i != retired.size(); ++i)
      m_retired_intervals.push_back(RetiredInterval(retired.interval(i), n));

    return n;
  }

  ARG::recomb_node_pair_t ARG::recombination(double time, Node *child,
    double cross_over_point)
    throw(null_event, null_child,
    interval_out_of_range,empty_interval)
  {
    if (child == 0) throw null_child();

    if (cross_over_point <= child->intervals().first_point())
      throw null_event();
    if (child->intervals().last_point() <= cross_over_point)
      throw null_event();

    Intervals left  = child->intervals().copy(0.0,cross_over_point);
    if (!m_keep_empty) left  = filter_contains_marker(left, m_conf);
    Intervals right = child->intervals().copy(cross_over_point,1.0);
    if (!m_keep_empty) right = filter_contains_marker(right, m_conf);

#if 0
    std::cout << "recombination -- child: " << child->intervals() << std::endl;
    std::cout << "recombination -- left: " << left << std::endl;
    std::cout << "recombination -- right: " << right << std::endl;
#endif

    RecombinationNode *n1 = new RecombinationNode(m_conf,time,child,left,
      cross_over_point, true);
    RecombinationNode *n2 = new RecombinationNode(m_conf,time,child,right,
      cross_over_point, false);
    m_node_pool.push_back(n1); m_node_pool.push_back(n2);

    return std::make_pair(n1,n2);
  }

  ARG::gene_conv_node_pair_t ARG::gene_conversion(double time, Node *child,
    double conversion_start,
    double conversion_end)
    throw(null_event, null_child,
    interval_out_of_range,empty_interval)
  {
    if (child == 0) throw null_child();

    if (conversion_start == conversion_end)
      throw null_event();

    if (conversion_end <= child->intervals().first_point())
      throw null_event();
    if (child->intervals().last_point() <= conversion_start)
      throw null_event();

    Intervals left  = child->intervals().copy(0.0, conversion_start)
      + child->intervals().copy(conversion_end, 1.0);
    if (!m_keep_empty) left  = filter_contains_marker(left, m_conf);

    Intervals right =
      child->intervals().copy(conversion_start, conversion_end);
    if (!m_keep_empty) right = filter_contains_marker(right, m_conf);

#if 0
    std::cout << "gene-conversion -- left: " << left << std::endl;
    std::cout << "gene-conversion -- right: " << right << std::endl;
#endif

    GeneConversionNode *n1 = new GeneConversionNode(m_conf,time,child,left,
      conversion_start,
      conversion_end,
      false);
    GeneConversionNode *n2 = new GeneConversionNode(m_conf,time,child,right,
      conversion_start,
      conversion_end,
      true);
    m_node_pool.push_back(n1); m_node_pool.push_back(n2);

    return std::make_pair(n1,n2);

  }

  namespace
  {
    using std::binary_function;
    struct starts_before :
    public binary_function<const Interval&,const Interval&,bool>
    {
      bool operator () (const Interval &i1, const Interval &i2) const
        { return i1.start() < i2.start(); }
    };
  };

  void ARG::sort_retired_intervals()
  {
    sort(m_retired_intervals.begin(), m_retired_intervals.end(),
      starts_before());
  }

  namespace
  {
    class state_printer : public std::unary_function<void,const Node*>
    {
      public:
        state_printer(std::ostream &os) : m_os(os) {};
        void operator () (const Node *n)
        {
          for (unsigned int s = 0; s < n->no_states(); ++s)
            m_os << n->state(s) << ' ';
          m_os << '\n';
        }

      private:
        std::ostream &m_os;
    };
  }

  void ARG::to_text(std::ostream &os) const
  {
    // header...
    os << "# markers: ";
    for (int i = 0; i < m_conf.no_markers(); ++i)
      os << m_conf.marker(i) << ' ' << m_conf.position(i) << ' ';
    os << '\n';

    // body...
    for_each(m_leaf_pool.begin(), m_leaf_pool.end(), state_printer(os));
  }

}
