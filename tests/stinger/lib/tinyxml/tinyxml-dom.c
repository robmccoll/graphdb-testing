// Copyright 2010 Christoph Gärtner
// Distributed under the Boost Software License, Version 1.0
// See <license.txt> or <http://www.boost.org/LICENSE_1_0.txt> for details

#include "tinyxml.h"
#include "tinyxml-dom.h"

dom_node_t dom_getElementById(
	dom_node_t root, const char *id)
{
	dom_node_t node = tx_find(root, NULL, TX_ATTRIBUTE, "id", id, 1);
	return node ? node->parent : NULL;
}

dom_node_t dom_nextSibling(
	dom_node_t node)
{
	return tx_next(node, NULL, 1, TX_ELEMENT | TX_TEXT);
}

dom_node_t dom_prevSibling(
	dom_node_t node)
{
	return tx_prev(node, NULL, 1, TX_ELEMENT | TX_TEXT);
}

dom_node_t dom_firstChild(
	dom_node_t node)
{
	return tx_next(node, node, 1, TX_ELEMENT | TX_TEXT);
}

const char *dom_getAttribute(
	dom_node_t node, const char *name)
{
	dom_node_t attribute = tx_find(node, NULL, TX_ATTRIBUTE, name, NULL, 0);
	return attribute ? attribute->value : NULL;
}

size_t dom_getElementsByTagName(
	dom_node_t root, const char *name, size_t max_nodes, dom_node_t *nodes)
{
	size_t count = 0;

	dom_node_t current = NULL;
	for(; count < max_nodes; ++count)
	{
		current = tx_find(root, current, TX_ELEMENT, name, NULL, 1);
		if(current) nodes[count] = current;
		else break;
	}

	return count;
}

_Bool dom_hasChildNodes(
	dom_node_t node)
{
	return tx_find(node, NULL, TX_ELEMENT | TX_TEXT, NULL, NULL, 0);
}

size_t dom_childNodes(
	dom_node_t node, size_t max_nodes, dom_node_t *nodes)
{
	size_t count = 0;

	dom_node_t current = NULL;
	for(; count < max_nodes; ++count)
	{
		current = tx_next(current, node, 1, TX_ELEMENT | TX_TEXT);
		if(current) nodes[count] = current;
		else break;
	}

	return count;
}
