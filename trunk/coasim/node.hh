/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__NODE_HH_INCLUDED
#define CORE__NODE_HH_INCLUDED

#ifndef CORE__CONFIGURATION_HH_INCLUDED
# include "configuration.hh"
#endif
#ifndef CORE__INTERVAL_HH_INCLUDED
# include "interval.hh"
#endif
#ifndef CORE__RETIRED_INTERVAL_HH_INCLUDED
# include "retired_interval.hh"
#endif

#ifndef STDEXCEPT_INCLUDED
# include <stdexcept>
# define STDEXCEPT_INCLUDED
#endif
#ifndef VECTOR_INCLUDED
# include <vector>
# define VECTOR_INCLUDED
#endif
#ifndef IOSTREAM_INCLUDED
# include <iostream>
# define IOSTREAM_INCLUDED
#endif
#ifndef VALARRAY_INCLUDED
# include <valarray>
# define VALARRAY_INCLUDED
#endif

namespace core
{

  class Marker;
  class Mutator;

  // -- Abstract class for ARG nodes --------------------------------------
  class Node
  {
    // explicitly remove chance of copying
    Node(const Node&);
    Node &operator = (const Node&);

    public:
      Node(const Configuration &conf, double time)
        : i_time(time), i_states(-1,conf.no_markers())
        {}
      Node(const Configuration &conf, double time, const Intervals &i)
        : i_time(time), i_intervals(i),  i_states(-1,conf.no_markers())
        {}
      virtual ~Node() {};

      double time() const { return i_time; }

      // the sub-intervals of the range [0,1) that connects this node to
      // a leaf node in the ARG -- if a point is not in one of these
      // intervals it is lost somewhere from here to the leaves.  As an
      // invariant, two points on the same interval correspond to the
      // same binary tree of the ARG.
      const Intervals &intervals() const { return i_intervals; }

      // calculate the "surface" (i.e. the total edge-length) of the
      // tree in point, with this node as root (zero if the point is not
      // in this node's intervals).
      virtual double surface_at_point(double point) const
        throw(std::out_of_range) = 0;

      virtual bool contains_point(double point) const
      {
        return i_intervals.contains_point(point);
      }

      // Prints the local tree to a stream
      void print_tree_at_point(std::ostream &os, double point) const
        throw(std::out_of_range)
      {
        print_tree_at_point(os, point, 0.0,false);
        os << ';';
      }
      // Prints the local tree to a stream
      virtual void print_tree_at_point(std::ostream &os, double point,
        double edge_length,
        bool print_edge) const
        throw(std::out_of_range) = 0;

      // Calculate the number of leaves hit by the binary tree in point
      // with root in this node.
      unsigned int leaves_at_point(double point) const throw(std::out_of_range)
        { return intervals().leaf_contacts(point); }

      // FIXME: I am not sure this is the access-protection for these...
      void initialize_marker(unsigned int idx, const Marker &m);
      virtual void mutate_marker(unsigned int idx, Mutator &m) = 0;

      unsigned int no_states()  const { return i_states.size(); }
      int state(unsigned int s) const throw(std::out_of_range)
      {
        if (i_states.size() <= s) throw std::out_of_range("s out of range");
        return i_states[s];
      }

    protected:
      void set_state(unsigned int s, int val)
        { i_states[s] = val; }

      // hack to work around C++'s crappy "don't access protected through
      // other than this" protection...
      static void set_state(Node *n, unsigned int s, int val)
        { n->set_state(s,val); }

    private:
      friend class ARG;

      double    i_time;
      Intervals i_intervals;
      std::valarray<int> i_states;
  };

  class LeafNode;
  class CoalescentNode;
  class RecombinationNode;
  class GeneConversionNode;

  // Exception thrown if a node is created with a 0-child
  struct null_child : public std::logic_error
  {
    null_child() : std::logic_error("null child.") {}
  };

  // Exception thrown if a recombination or gene conversion
  // falls outside an interval, and the even should just be ignored
  struct null_event : public std::exception {};

  class ARG
  {
    public:

      // -- Initialization and book-keeping -----------------------------------
      ARG(const Configuration &conf, bool keep_empty_intervals = false)
        : i_conf(conf), i_keep_empty(keep_empty_intervals),
        i_no_leaves(0)
        {}

      // Cleanup.  Destroying the ARG also deletes all nodes in it, so
      // don't keep any pointers to them around after this deletion.
      ~ARG();

      // -- Factory methods for building the ARG ------------------------------
      LeafNode *leaf() throw();
      CoalescentNode *coalescence(double time, Node *left, Node *right)
        throw(null_child);

      // these methods return a pair of new nodes, or throws a null_event
      // exception if the event should be ignored.

      typedef std::pair<RecombinationNode*,RecombinationNode*>
        recomb_node_pair_t;
      recomb_node_pair_t recombination(double time, Node *child,
        double cross_over_point)
        throw(null_event, null_child,
        interval_out_of_range,empty_interval);
      typedef std::pair<GeneConversionNode*,GeneConversionNode*>
        gene_conv_node_pair_t;
      gene_conv_node_pair_t gene_conversion(double time, Node *child,
        double conversion_start,
        double conversion_end)
        throw(null_event, null_child,
        interval_out_of_range,empty_interval);

      const std::vector<RetiredInterval> &retired_intervals() const
        { return i_retired_intervals; }
      void sort_retired_intervals();

      unsigned int no_nodes() const
        { return i_leaf_pool.size() + i_node_pool.size(); }

      void to_text(std::ostream &os) const;

      const std::vector<Node*> &inner_nodes() const { return i_node_pool; }
      const std::vector<Node*> &leaves()      const { return i_leaf_pool; }

    private:
      // disable these
      ARG(const ARG&);
      ARG &operator = (const ARG&);

      const Configuration &i_conf;
      bool i_keep_empty;
      size_t i_no_leaves;

      // pools of nodes -- FIXME: can be handled more efficiently...
      std::vector<Node*> i_leaf_pool;
      std::vector<Node*> i_node_pool;

      std::vector<RetiredInterval> i_retired_intervals;
  };

  class LeafNode : public Node
  {
    friend LeafNode *ARG::leaf();
    LeafNode(const Configuration &conf, unsigned int id)
      : Node(conf,0.0), i_id(id)
      {}

    virtual double surface_at_point(double point) const
      throw(std::out_of_range);
    virtual void print_tree_at_point(std::ostream &os, double point,
      double edge_length,
      bool print_edge) const
      throw(std::out_of_range);
    virtual void mutate_marker(unsigned int idx, Mutator &m);

    unsigned int i_id;
  };

  // WARNING: None of the following classes checks whether they are
  // initialized with a null-child, but will crash if that is the
  // case.  They should only be created with the corresponding factory
  // method anyway, and it checks for it, so *DON'T* create these
  // objects in any other way!

  class CoalescentNode : public Node
  {
    public:
      const Node * const left_child()  const { return i_left; }
      const Node * const right_child() const { return i_right; }

      const Intervals &retired_intervals() const
      {
        return i_retired_intervals;
      }

    private:
      friend CoalescentNode *ARG::coalescence(double,Node*,Node*);
      CoalescentNode(const Configuration &conf, double time,
        Node *left, Node *right,
        const Intervals &is, const Intervals &ris)
        : Node(conf,time,is), i_left(left), i_right(right),
        i_retired_intervals(ris), i_conf(conf)
        {}

      virtual double surface_at_point(double point) const
        throw(std::out_of_range);

      virtual bool contains_point(double point) const
      {
        return Node::contains_point(point) or
          i_retired_intervals.contains_point(point);
      }

      virtual void print_tree_at_point(std::ostream &os, double point,
        double edge_length,
        bool print_edge) const
        throw(std::out_of_range);
      virtual void mutate_marker(unsigned int idx, Mutator &m);

      Node *const i_left;
      Node *const i_right;
      Intervals i_retired_intervals;
      const Configuration &i_conf;

  };

  class RecombinationNode : public Node
  {
    public:
      double cross_over_point()  const { return i_cross_over_point; }
      const Node * const child() const { return i_child; }

    private:
      friend ARG::recomb_node_pair_t ARG::recombination(double,Node*,double);
      RecombinationNode(const Configuration &conf,
        double time, Node *child, const Intervals &is,
        double cross_over_point, bool is_left)
        : Node(conf,time,is), i_child(child),
        i_cross_over_point(cross_over_point), i_is_left(is_left)
        {}

      virtual double surface_at_point(double point) const
        throw(std::out_of_range);
      virtual void print_tree_at_point(std::ostream &os, double point,
        double edge_length,
        bool print_edge) const
        throw(std::out_of_range);
      virtual void mutate_marker(unsigned int idx, Mutator &m);

      Node *const i_child;
      double i_cross_over_point;
      bool i_is_left;

  };

  class GeneConversionNode : public Node
  {
    public:
      double conversion_start() const { return i_conversion_start; }
      double conversion_end()   const { return i_conversion_end; }
      const Node * const child() const { return i_child; }

    private:
      friend ARG::gene_conv_node_pair_t ARG::gene_conversion(double,Node*,double,double);
      GeneConversionNode(const Configuration &conf,
        double time, Node *child, const Intervals &is,
        double conversion_start, double conversion_end,
        bool is_inside)
        : Node(conf,time,is), i_child(child),
        i_conversion_start(conversion_start), i_conversion_end(conversion_end),
        i_is_inside(is_inside)
        {}

      virtual double surface_at_point(double point) const
        throw(std::out_of_range);
      virtual void print_tree_at_point(std::ostream &os, double point,
        double edge_length,
        bool print_edge) const
        throw(std::out_of_range);
      virtual void mutate_marker(unsigned int idx, Mutator &m);

      Node *const i_child;
      double i_conversion_start, i_conversion_end;
      bool i_is_inside;
  };

}
#endif
