Working with modified peptides
------------------------------

When trying to score the localization of a modification on a peptide's sequence,
it is often necessary to enumerate all possible localizations of the modification
and calculate a score. pyAscore has an internal class for doing this enumeration,
and we make that interface available to users to allow convenient iteration
over the theoretical fragments of modified peptides.

Iteration over all fragments
############################

Basic iteration requires a peptide sequence (e.g. 'ASMTK'), the mass of a modification
of interest (e.g. 79.9663), the amino acids that the modification can fall on (e.g STY),
and the number of modifications of interest. All this can be initialized for a given
peptide with the following lines.

.. code-block:: python

   pep = PyModifiedPeptide("STY", 79.9663)
   pep.consume_peptide("ASMTK", 1)

The **PyModifiedPeptide** object is where we will get our graphs which allow iteration
over the theoretical fragments from permutations of modified amino acids. In this case,
there are 2 permutations, AS[79.9663]MTK and ASMT[79.9663]K. If we would like to iterate
over the b fragments of charge 1, we can do that by generating the b+ fragment graph
and iterating using the **iter_permutations** and **iter_fragments** methods.

.. code-block:: python

   b_graph = pep.get_fragment_graph("b", 1)
   for perm in b_graph.iter_permutations():
       print(perm.get_signature())
       for mz, label in graph.iter_fragments():
           print(mz, "<-", label)

.. code-block:: none

   # Output:
   [1, 0]
   72.0449 <- b1+
   239.043 <- b2+
   370.084 <- b3+
   471.131 <- b4+
   [0, 1]
   72.0449 <- b1+
   159.077 <- b2+
   290.117 <- b3+
   471.131 <- b4+

Notice that the full peptide mass is currently not included. This is because this "fragment"
can't be used for localization.

When recording scores, we tend to build a dictionary that uses the signature as a key.
This allows us to track additive scores, such as the counts pyAscore uses, and iterate
over graphs independently of each other. Note that one of the big reasons we do that
is that iteration over b type permutations is in the opposite direction to the y type
permutations.

.. code-block:: python

   print("b type iteration:")
   b_graph = pep.get_fragment_graph("b", 1)
   for perm in b_graph.iter_permutations():
       print(perm.get_signature())

   print("y type iteration:")
   y_graph = pep.get_fragment_graph("y", 1)
   for perm in y_graph.iter_permutations():
       print(perm.get_signature())

.. code-block:: none

   # Output
   b type iteration:
   [1, 0]
   [0, 1]
   y type iteration:
   [0, 1]
   [1, 0]

.. autoclass:: pyascore.PyModifiedPeptide
   :members:

.. autoclass:: pyascore.PyFragmentGraph
   :members:
