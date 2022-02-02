.. autoclass:: {{ objname }}
   
    {% block methods %}
    {% if methods %}
    .. rubric:: Methods

    .. autosummary::
        :toctree: 
        
        {% for item in methods %}
            {%- if item != '__init__' %}
            ~{{ fullname }}.{{ item }}
            {%- endif -%}
        {%- endfor %}
        {% endif %}
        {% endblock %}
   
