// Copyright 2010 Christoph Gärtner
// Distributed under the Boost Software License, Version 1.0
// See <license.txt> or <http://www.boost.org/LICENSE_1_0.txt> for details

#include "tinyxml.h"

#include <ctype.h>
#include <string.h>

struct tx_node *tx_next(
	struct tx_node *node, struct tx_node *ancestor, _Bool child,
	enum tx_types type)
{
	if(!ancestor) ancestor = node->parent;

	for(++node; node->type && node->parent >= ancestor; ++node)
	{
		if((node->type & type) && (!child || node->parent == ancestor))
			return node;
	}

	return NULL;
}

struct tx_node *tx_prev(
	struct tx_node *node, struct tx_node *ancestor, _Bool child,
	enum tx_types type)
{
	if(!ancestor) ancestor = node->parent;

	for(--node; node > ancestor; --node)
	{
		if((node->type & type) && (!child || node->parent == ancestor))
			return node;
	}

	return NULL;
}

struct tx_node *tx_find(
	struct tx_node *root, struct tx_node *current, enum tx_types type,
	const char *name, const char *value, _Bool deep)
{
	if(!current) current = root;

	for(;;)
	{
		current = tx_next(current, root, !deep, type);
		if(!current || (
			(!name || strcmp(current->name, name) == 0) &&
			(!value || strcmp(current->value, value) == 0)))
			return current;
	}
}

struct tx_node *tx_get(
	struct tx_node *root, struct tx_node *current, size_t path_len,
	const char **path)
{
	enum tx_types type = (path[path_len - 1][0] == '@') ? TX_ATTRIBUTE :
		(strcmp(path[path_len - 1], "#text") == 0) ? TX_TEXT : TX_ELEMENT;

	for(;;)
	{
		current = tx_find(root, current, type,
			path[path_len - 1] + (type == TX_ATTRIBUTE ? 1 : 0), NULL, 1);

		if(!current || path_len == 1)
			return current;

		struct tx_node *ancestor = current->parent;
		size_t pos = path_len - 2;

		for(; ancestor > root && strcmp(path[pos], ancestor->name) == 0;
			ancestor = ancestor->parent, --pos)
		{
			if(pos == 0)
				return current;
		}
	}
}

#define call(LABEL) \
	do { ++data; goto LABEL; } while(0)

#define put_node(...) \
	do { if(!count--) return data; put_node_(&nodes, __VA_ARGS__); } while(0)

static inline void put_node_(
	struct tx_node **nodes, enum tx_types type, struct tx_node *parent,
	const char *name, const char *value)
{
	(*nodes)->type = type;
	(*nodes)->parent = parent;
	(*nodes)->name = name;
	(*nodes)->value = value;
	++*nodes;
}

char *tx_parse(
	char *data, size_t count, struct tx_node *nodes)
{
	char *marks[3] = { data };
	struct tx_node *parent = nodes;
	put_node(TX_EOF, NULL, NULL, NULL);

	parse_text:
	{
		while(*data && *data != '<') ++data;

		char c = *data;

		if(data > marks[0])
		{
			put_node(TX_TEXT, parent, "#text", marks[0]);
			*data = 0;

			if(!parent->value)
				parent->value = marks[0];
		}

		if(c == '<') call(parse_element);

		put_node(TX_EOF, NULL, NULL, NULL);
		return NULL;
	}

	parse_element:
	{
		if(*data == '/') call(parse_closing_element);

		while(isspace(*data)) ++data;

		if(isalnum(*data))
		{
			marks[0] = data;
			call(parse_element_name);
		}

		return data;
	}

	parse_closing_element:
	{
		while(isspace(*data)) ++data;

		char *start = data;
		while(isalnum(*data)) ++data;
		while(isspace(*data)) ++data;

		if(*data != '>')
			return data;

		*data = 0;

		if(!parent->name || strcmp(parent->name, start))
			return start;

		parent = parent->parent;

		marks[0] = data + 1;
		call(parse_text);
	}

	parse_element_name:
	{
		while(isalnum(*data)) ++data;

		put_node(TX_ELEMENT, parent, marks[0], NULL);
		parent = nodes - 1;

		char c = *data;
		*data = 0;

		if(c == '>')
		{
			marks[0] = data + 1;
			call(parse_text);
		}
		else if(c == '/')
		{
			if(*++data == '>')
			{
				parent = parent->parent;
				marks[0] = data + 1;
				call(parse_text);
			}

			return data;
		}
		else if(isspace(c))
		{
			call(parse_attributes);
		}

		return data;
	}

	parse_attributes:
	{
		while(isspace(*data)) ++data;

		if(isalnum(*data))
		{
			marks[0] = data;
			call(parse_attribute_name);
		}
		else if(*data == '>')
		{
			marks[0] = data + 1;
			call(parse_text);
		}
		else if(*data == '/')
		{
			if(*++data == '>')
			{
				parent = parent->parent;
				marks[0] = data + 1;
				call(parse_text);
			}

			return data;
		}

		return data;
	}

	parse_attribute_name:
	{
		while(isalnum(*data)) ++data;

		marks[1] = data;

		while(isspace(*data)) ++data;

		if(*data == '=') call(parse_attribute_value);

		return data;
	}


	parse_attribute_value:
	{
		while(isspace(*data)) ++data;

		marks[2] = data + 1;

		char c = *data;

		if(c != '"' && c != '\'')
			return data;

		++data;

		for(; *data; ++data)
		{
			if(*data == '\\')
			{
				if(*++data == 0)
					return data;
			}
			else if(*data == c)
			{
				put_node(TX_ATTRIBUTE, parent, marks[0], marks[2]);

				*marks[1] = 0;
				*data = 0;

				call(parse_attributes);
			}
		}

		return data;
	}
}
