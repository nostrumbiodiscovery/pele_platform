
The reStructuredText Cheat Sheet: Basic Syntax
==================================================


Headings
**********
When writing a title, make sure to write as many symbols (=,-,~,*,^) as there are characters in the title.
Titles are underlined (but also can be overlined) with one of the previous symbols. Any non-alphanumeric character can be used but the usual convention is:

* ``#`` with overline, for parts.
* ``*`` with overline, for chapters.
* ``=`` for sections.
* ``-`` for subsections.
* ``^`` for subsubsections.
* ``"`` for paragraphs.

Sections
**********
::

    ================================================================================
    Document title
    ================================================================================

    First level
    -----------

    Second level
    ++++++++++++

    Third level
    ***********

    Fourth level
    ~~~~~~~~~~~~

Inline Markup
***************

+----------------------------------+--------------------------------+
| ``*italic*``                     | *italic*                       |
+----------------------------------+--------------------------------+
| ``**bold**``                     | **bold**                       |
+----------------------------------+--------------------------------+
| ``reference_``                   | reference_                     |
+----------------------------------+--------------------------------+
|  ``footnote reference [1]_``     | footnote reference [1]_        |
+----------------------------------+--------------------------------+


List and bullets
********************
+---------------------------------------+---------------------------------+
|  ``* This is a bulleted list``        | * This is a bulleted list       |
|                                       |                                 |
|  ``* Second item``                    | * Second item                   |
+---------------------------------------+---------------------------------+
|  ``1. This is a numbered list``       | 1. This is a numbered list      |
|                                       |                                 | 
|  ``2. Second item``                   | 2. Second item                  |
+---------------------------------------+---------------------------------+
|``#. This is also a numbered list.``   | #. This is also a numbered list |
|                                       |                                 |
|``#. Second item``                     | #. Second item                  |
+---------------------------------------+---------------------------------+ 


