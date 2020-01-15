---
layout: page
title: Syllabus
img: lab.png # Add image post (optional)
permalink: syllabus
sidebar: true
---

<table>
<tr>
    <th><b>Date</b></th>
    <th><b>Week</b></th>
    <th><b>Topic</b></th>
    <th><b>Reading</b></th>
    <th><b>Tasks Due</b></th>
</tr>
{% for day in site.data.syllabus %}
<tr>
    <td>{{day.date}}</td>
    <td>{{day.week}}</td>
    <td>{{day.topic}}</td>
    <td>{{day.reading}}</td>
    {% if day.due %}
    <td>{{day.due}}</td>
    {% else %}
    <td> -- </td>
    {% endif %}
    {% if day.notes %}
    <td><a href="https://rpdata.caltech.edu/courses/bige105/protected/{{day.notes}}">
    PDF </a></td>
    {% else %}
    <td> -- </td>
    {% endif %}
</tr>
{% endfor %}
</table>
