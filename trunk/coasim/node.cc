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

using namespace core;

void
core::Node::initialize_marker(unsigned int idx, const Marker &m)
{
  if (i_states.size() <= idx)
    throw std::out_of_range("marker index out of range");
  i_states[idx] = m.default_value();
}


double
core::LeafNode::surface_at_point(double point) const
throw(std::out_of_range)
{
  if (point < 0 or 1.0 <= point)
    throw std::out_of_range("Point out of range [0,1).");
  return 0.0;
}


void
core::LeafNode::print_tree_at_point(std::ostream &os, double point,
double edge_length,
bool print_edge) const
throw(std::out_of_range)
{
  if (point < 0 or 1.0 <= point)
    throw std::out_of_range("Point out of range [0,1).");
  os << '\'' << i_id << '\'';
  if (print_edge) os << " : " << edge_length;
}


void
core::LeafNode::mutate_marker(unsigned int idx, Mutator &m)
{
  // no mutations out of leaf
}


double
core::CoalescentNode::surface_at_point(double point) const
throw(std::out_of_range)
{
  // NB! don't check if this node contains it -- it could be
  // retired, if that is the case it's children will contain it.
  // if this node does not contain it, neither of it's children
  // will, so it still works out (althoug we call recursively at
  // little more than strictly necessary)
  double surface = 0.0;
  if (i_left->intervals().contains_point(point))
  {
    surface += i_left->surface_at_point(point);
    surface += time() - i_left->time();
  }
  if (i_right->intervals().contains_point(point))
  {
    surface += i_right->surface_at_point(point);
    surface += time() - i_right->time();
  }
  return surface;
}


void
core::CoalescentNode::print_tree_at_point(std::ostream &os, double point,
double edge_length,
bool print_edge) const
throw(std::out_of_range)
{
  double left_dist  = time() - i_left->time();
  double right_dist = time() - i_right->time();

  if (i_left->intervals().contains_point(point)
    and i_right->intervals().contains_point(point))
  {
    os << '(';
    i_left->print_tree_at_point(os, point, left_dist, true);
    os << ',';
    i_right->print_tree_at_point(os, point, right_dist, true);
    os << ')';
    if (print_edge) os << " : " << edge_length;

  }
  else
  {
    if (i_left->intervals().contains_point(point))
      i_left->print_tree_at_point(os, point,
        edge_length+left_dist,
        print_edge);
    if (i_right->intervals().contains_point(point))
      i_right->print_tree_at_point(os, point,
        edge_length+right_dist,
        print_edge);
  }
}


void
core::CoalescentNode::mutate_marker(unsigned int idx, Mutator &m)
{
  if (! (idx < no_states()) )
    throw std::out_of_range("marker index out of range");

  double point = i_conf.position(idx);

  if (i_left->intervals().contains_point(point))
  {
    set_state(i_left, idx, m.mutate(*this,*i_left,state(idx)));
    i_left->mutate_marker(idx,m);
  }

  if (i_right->intervals().contains_point(point))
  {
    set_state(i_right, idx, m.mutate(*this,*i_right,state(idx)));
    i_right->mutate_marker(idx,m);
  }
}


double
core::RecombinationNode::surface_at_point(double point) const
throw(std::out_of_range)
{
  // no need to check here, it is the parents responsibility to
  //check that
  //if (! intervals().contains_point(point)) return 0.0;

  double surface = 0.0;
  if (i_child->intervals().contains_point(point))
  {
    surface += i_child->surface_at_point(point);
    surface += time() - i_child->time();
  }
  return surface;
}


void
core::RecombinationNode::print_tree_at_point(std::ostream &os, double point,
double edge_length,
bool print_edge) const
throw(std::out_of_range)
{
  double d = time() - i_child->time();
  i_child->print_tree_at_point(os, point, edge_length+d, print_edge);
}


void
core::RecombinationNode::mutate_marker(unsigned int idx, Mutator &m)
{
  if (! (idx < no_states()) )
    throw std::out_of_range("marker index out of range");

  set_state(i_child, idx, m.mutate(*this,*i_child,state(idx)));
  i_child->mutate_marker(idx,m);
}


double
core::GeneConversionNode::surface_at_point(double point) const
throw(std::out_of_range)
{
  // no need to check here, it is the parents responsibility to
  //check that
  //if (! intervals().contains_point(point)) return 0.0;

  double surface = 0.0;
  if (i_child->intervals().contains_point(point))
  {
    surface += i_child->surface_at_point(point);
    surface += time() - i_child->time();
  }
  return surface;
}


void
core::GeneConversionNode::print_tree_at_point(std::ostream &os, double point,
double edge_length,
bool print_edge) const
throw(std::out_of_range)
{
  double d = time() - i_child->time();
  i_child->print_tree_at_point(os, point, edge_length+d, print_edge);
}


void
core::GeneConversionNode::mutate_marker(unsigned int idx, Mutator &m)
{
  if (! (idx < no_states()) )
    throw std::out_of_range("marker index out of range");

  set_state(i_child, idx, m.mutate(*this,*i_child,state(idx)));
  i_child->mutate_marker(idx,m);
}


ARG::~ARG()
{
  std::vector<Node*>::iterator itr;
  for (itr = i_leaf_pool.begin(); itr != i_leaf_pool.end(); ++itr)
    delete *itr;
  for (itr = i_node_pool.begin(); itr != i_node_pool.end(); ++itr)
    delete *itr;
}


LeafNode *ARG::leaf() throw()
{
  LeafNode *n = new LeafNode(i_conf, i_no_leaves);
  n->i_intervals.add(0.0,1.0,1);
  i_leaf_pool.push_back(n);
  ++i_no_leaves;

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
  if (!i_keep_empty) merged = filter_contains_marker(merged, i_conf);

#if 0
  std::cout << "coalescence -- left: " << left->intervals() << std::endl;
  std::cout << "coalescence -- right: " << right->intervals() << std::endl;
  std::cout << "coalescence -- merged: " << merged << std::endl;
#endif

  for (int i = 0; i < merged.size(); ++i)
  {
    if (merged.interval(i).leaf_contacts() < i_no_leaves)
      non_retired.add(merged.interval(i));
    else if (merged.interval(i).leaf_contacts() == i_no_leaves)
      retired.add(merged.interval(i));
    else
      assert(false);
  }

  CoalescentNode *n = new CoalescentNode(i_conf,time,left,right,
    non_retired, retired);

  i_node_pool.push_back(n);

  for (int i = 0; i != retired.size(); ++i)
    i_retired_intervals.push_back(RetiredInterval(retired.interval(i), n));

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
  if (!i_keep_empty) left  = filter_contains_marker(left, i_conf);
  Intervals right = child->intervals().copy(cross_over_point,1.0);
  if (!i_keep_empty) right = filter_contains_marker(right, i_conf);

#if 0
  std::cout << "recombination -- child: " << child->intervals() << std::endl;
  std::cout << "recombination -- left: " << left << std::endl;
  std::cout << "recombination -- right: " << right << std::endl;
#endif

  RecombinationNode *n1 = new RecombinationNode(i_conf,time,child,left,
    cross_over_point, true);
  RecombinationNode *n2 = new RecombinationNode(i_conf,time,child,right,
    cross_over_point, false);
  i_node_pool.push_back(n1); i_node_pool.push_back(n2);

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
  if (!i_keep_empty) left  = filter_contains_marker(left, i_conf);

  Intervals right =
    child->intervals().copy(conversion_start, conversion_end);
  if (!i_keep_empty) right = filter_contains_marker(right, i_conf);

#if 0
  std::cout << "gene-conversion -- left: " << left << std::endl;
  std::cout << "gene-conversion -- right: " << right << std::endl;
#endif

  GeneConversionNode *n1 = new GeneConversionNode(i_conf,time,child,left,
    conversion_start,
    conversion_end,
    false);
  GeneConversionNode *n2 = new GeneConversionNode(i_conf,time,child,right,
    conversion_start,
    conversion_end,
    true);
  i_node_pool.push_back(n1); i_node_pool.push_back(n2);

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
  sort(i_retired_intervals.begin(), i_retired_intervals.end(),
    starts_before());
}


namespace
{
  class state_printer : public std::unary_function<void,const Node*>
  {
    public:
      state_printer(std::ostream &os) : i_os(os) {};
      void operator () (const Node *n)
      {
        for (unsigned int s = 0; s < n->no_states(); ++s)
          i_os << n->state(s) << ' ';
        i_os << '\n';
      }

    private:
      std::ostream &i_os;
  };
}


void ARG::to_text(std::ostream &os) const
{
  // header...
  os << "# markers: ";
  for (int i = 0; i < i_conf.no_markers(); ++i)
    os << i_conf.marker(i) << ' ' << i_conf.position(i) << ' ';
  os << '\n';

  // body...
  for_each(i_leaf_pool.begin(), i_leaf_pool.end(), state_printer(os));
}
