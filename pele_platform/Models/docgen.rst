YAML Model
==========

{% for field in fields %}
{{field.title}}
########################################################################

{% if field.description is none %}
Description missing.
{% else %}
{{field.description}}
{% endif %}

{% if field.from_simulation_params %}
:YAML key: None, value calculated in simulation params
{% else %}
{% if field.value_from is none %}
:YAML key: {{field.alias}}
{% else %}
:YAML key: None, value from `{{ field.value_from }}` parser
{% endif %}
{% endif %}
:Parser key: {{field.name}}
:Type: ``{{field.type_}}``
{% if field.default is not none %}
:Default value: ``{{ field.default }}``
{% endif %}
{% if field.tests_value is not none %}
:Tests value: ``{{ field.tests_value }}``
{% endif %}

{% if field.value_from_simulation_params is not none %}
If value is {% if field.can_be_falsy %}``None``{% else %} falsy{% endif %},
it will fall back to ``{{ field.value_from_simulation_params }}`` software setting.
{% if field.simulation_params_default is not none %}If there is no software setting,
it will default to ``{{ field.simulation_params_default }}``.
{% endif %}
{% endif %}

{% if field.validators %}
Processors:{% for validator in field.validators %}
 * {{validator}}
{% endfor %}
{% endif %}

{% endfor %}