Online help system
==================

Most of the help information contained in this document and *the simuPOP
reference manual* is available from command line. For example, after you install
and import the simuPOP module, you can use ``help(Population.addInfoField)``\ to
view the help information of member function ``addInfoField`` of class
:class:`Population`.

**Example**: *Getting help using the \texttt{help()} function*

::

   >>> import simuPOP as sim
   >>> help(sim.Population.addInfoFields)
   Help on built-in function Population_addInfoFields in module simuPOP._simuPOP_std:

   Population_addInfoFields(...)
       Usage:

           x.addInfoFields(fields, init=0)

       Details:

           Add a list of information fields fields to a population and
           initialize their values to init. If an information field alreay
           exists, it will be re-initialized.


   now exiting runScriptInteractively...

`Download help.py <help.py>`_

It is important that you understand that

* The constructor of a class is named ``__init__`` in Python. That is to say,
  you should use the following command to display the help information of the
  constructor of class :class:`Population`:

  ::

     >>> help(Population.__init__)

* Some classes are derived from other classes and have access to member
  functions of their base classes. For example, class :class:`Population` and
  :class:`Individual` are both derived from class :class:`GenoStruTrait`.
  Therefore, you can use all :class:`GenoStruTrait` member functions from these
  classes.

  In addition, the constructor of a derived class also calls the constructor of
  its base class so you may have to refer to the base class for some parameter
  definitions. For example, parameters ``begin, end, step, at``\ etc are shared by
  all operators, and are explained in details only in class ``BaseOperator.``


