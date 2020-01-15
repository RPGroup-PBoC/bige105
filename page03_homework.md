---
layout: page
title: Homework
img: code.png # Add image post (optional)
permalink: homework
sidebar: true
---

---


<table>
<tr>
<th> <b>Homework</b></th>
<th> <b>File</b></th>
<th> <b>Associated Reading</b></th>
<th> <b> Due Date</b> </th>
</tr>

{% for hwk in site.data.homework %}
<tr>
    <td>{{hwk.number}}</td>
    <td> <a href="http://www.rpgroup.caltech.edu/bige105/hwk/{{hwk.pset}}"> Problem Set </a></td>
    {% if hwk.reading %}
    <td> <a href="http://www.rpgroup.caltech.edu/bige105/hwk/{{hwk.reading}}"> Paper </a></td>
    {% else %}
    <td> -- </td>
    {% endif %}
    <td> {{hwk.due_date}} </td>
<tr>
{%endfor%}
</table>
