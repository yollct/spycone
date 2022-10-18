{% if functions %}

.. autosummary::
    :toctree: .
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}