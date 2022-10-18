{{ objname | escape | underline }}

.. rubric:: Description

.. currentmodule:: {{ objname }}


{% if functions %}
.. rubric:: Functions

.. autosummary::
    :toctree: .
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}